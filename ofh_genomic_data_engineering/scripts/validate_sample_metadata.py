from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

try:
    from scripts.metrics_to_parquet import semantic_table_sha256
except ModuleNotFoundError:
    from metrics_to_parquet import semantic_table_sha256


REQUIRED_METADATA_COLUMNS = ("sample", "pop", "super_pop", "gender")
CONTRACT = "sample-metadata-join-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("utf-8")
    return hashlib.sha1(header + data).hexdigest()


def ordered_ids_sha256(values: list[str]) -> str:
    raw = "\n".join(values).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def canonical_plink_iid(iid: str) -> str:
    parts = iid.split("_")
    if len(parts) != 2 or not parts[0] or parts[0] != parts[1]:
        raise ValueError(
            f"expected PLINK --double-id IID '<sample>_<sample>', observed {iid!r}"
        )
    return parts[0]


def read_psam(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=r"\s+", dtype=str, keep_default_na=False)
    frame.columns = [column.lstrip("#") for column in frame.columns]
    if "IID" not in frame.columns:
        raise ValueError("PSAM is missing IID")
    return frame


def read_metadata(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0)
    missing = [column for column in REQUIRED_METADATA_COLUMNS if column not in header.columns]
    if missing:
        raise ValueError(f"sample metadata is missing required columns: {missing}")
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=list(REQUIRED_METADATA_COLUMNS),
    )


def validate_sample_metadata(
    *,
    psam_path: Path,
    metadata_path: Path,
    expected_git_blob_sha1: str,
    source_url: str,
    output_parquet: Path,
    output_json: Path,
) -> dict[str, object]:
    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, detail: object) -> None:
        checks.append({"name": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    observed_blob = git_blob_sha1(metadata_path)
    record(
        "metadata_git_blob_sha1",
        observed_blob == expected_git_blob_sha1,
        {"expected": expected_git_blob_sha1, "observed": observed_blob},
    )

    try:
        metadata = read_metadata(metadata_path)
        metadata_read_error = None
    except ValueError as exc:
        metadata = pd.DataFrame(columns=REQUIRED_METADATA_COLUMNS)
        metadata_read_error = str(exc)
    record("metadata_required_columns", metadata_read_error is None, metadata_read_error or "present")

    try:
        psam = read_psam(psam_path)
        psam_read_error = None
    except ValueError as exc:
        psam = pd.DataFrame(columns=["IID"])
        psam_read_error = str(exc)
    record("psam_iid_column", psam_read_error is None, psam_read_error or "present")

    metadata_ids = metadata.get("sample", pd.Series(dtype="string")).astype(str)
    metadata_nonempty = bool(len(metadata_ids) > 0 and metadata_ids.str.len().gt(0).all())
    record("metadata_sample_ids_nonempty", metadata_nonempty, {"rows": len(metadata_ids)})
    record(
        "metadata_sample_ids_unique",
        bool(metadata_ids.is_unique),
        {"rows": len(metadata_ids), "unique": int(metadata_ids.nunique())},
    )

    plink_iids = psam.get("IID", pd.Series(dtype="string")).astype(str)
    record(
        "plink_iids_unique",
        bool(len(plink_iids) > 0 and plink_iids.is_unique),
        {"rows": len(plink_iids), "unique": int(plink_iids.nunique())},
    )

    canonical_ids: list[str] = []
    malformed: list[str] = []
    for iid in plink_iids.tolist():
        try:
            canonical_ids.append(canonical_plink_iid(iid))
        except ValueError:
            malformed.append(iid)
    record(
        "plink_iid_canonical_form",
        not malformed and len(canonical_ids) == len(plink_iids),
        {"malformed": malformed[:20], "malformed_count": len(malformed)},
    )
    record(
        "canonical_sample_ids_unique",
        bool(canonical_ids and len(set(canonical_ids)) == len(canonical_ids)),
        {"rows": len(canonical_ids), "unique": len(set(canonical_ids))},
    )

    metadata_lookup = metadata.set_index("sample", drop=False) if metadata_ids.is_unique else None
    missing_samples: list[str] = []
    if metadata_lookup is not None and len(canonical_ids) == len(plink_iids):
        missing_samples = [sample for sample in canonical_ids if sample not in metadata_lookup.index]
    else:
        missing_samples = canonical_ids.copy()
    matched_count = len(canonical_ids) - len(missing_samples)
    coverage = matched_count / len(canonical_ids) if canonical_ids else 0.0
    record(
        "metadata_full_coverage",
        bool(canonical_ids and not missing_samples and coverage == 1.0),
        {
            "plink_samples": len(canonical_ids),
            "matched": matched_count,
            "missing": missing_samples[:20],
            "missing_count": len(missing_samples),
            "coverage": coverage,
        },
    )

    joined: pd.DataFrame | None = None
    if metadata_lookup is not None and canonical_ids and not missing_samples and not malformed:
        rows = []
        for iid, sample_id in zip(plink_iids.tolist(), canonical_ids, strict=True):
            meta = metadata_lookup.loc[sample_id]
            rows.append(
                {
                    "IID": iid,
                    "sample": sample_id,
                    "pop": str(meta["pop"]),
                    "super_pop": str(meta["super_pop"]),
                    "gender": str(meta["gender"]),
                }
            )
        joined = pd.DataFrame(rows, columns=["IID", "sample", "pop", "super_pop", "gender"])

    joined_fields_complete = bool(
        joined is not None
        and not joined.empty
        and joined[["IID", "sample", "pop", "super_pop", "gender"]]
        .astype(str)
        .apply(lambda column: column.str.len().gt(0).all())
        .all()
    )
    record(
        "joined_fields_nonempty",
        joined_fields_complete,
        {"rows": 0 if joined is None else len(joined)},
    )
    ordered_identity = bool(
        joined is not None
        and joined["IID"].astype(str).tolist() == plink_iids.astype(str).tolist()
        and joined["sample"].astype(str).tolist() == canonical_ids
    )
    record("joined_order_identity", ordered_identity, {"rows": 0 if joined is None else len(joined)})

    passed = all(check["status"] == "PASS" for check in checks)
    semantic_hash = None
    if passed and joined is not None:
        output_parquet.parent.mkdir(parents=True, exist_ok=True)
        joined.to_parquet(output_parquet, index=False, compression="zstd")
        semantic_hash = semantic_table_sha256(joined)

    payload = {
        "status": "PASS" if passed else "FAIL",
        "contract": CONTRACT,
        "source": {
            "url": source_url,
            "git_blob_sha1": observed_blob,
            "sha256": sha256(metadata_path),
            "row_count": int(len(metadata)),
        },
        "join": {
            "plink_sample_count": int(len(plink_iids)),
            "matched_sample_count": int(matched_count),
            "coverage": float(coverage),
            "canonical_sample_ids_sha256": ordered_ids_sha256(canonical_ids) if canonical_ids else None,
            "output_semantic_sha256": semantic_hash,
            "population_counts": joined["pop"].value_counts().sort_index().to_dict() if joined is not None else {},
            "super_population_counts": joined["super_pop"].value_counts().sort_index().to_dict() if joined is not None else {},
            "gender_counts": joined["gender"].value_counts().sort_index().to_dict() if joined is not None else {},
        },
        "checks": checks,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--psam", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--expected-git-blob-sha1", required=True)
    p.add_argument("--source-url", required=True)
    p.add_argument("--output-parquet", required=True)
    p.add_argument("--output-json", required=True)
    return p


def main() -> int:
    args = parser().parse_args()
    payload = validate_sample_metadata(
        psam_path=Path(args.psam),
        metadata_path=Path(args.metadata),
        expected_git_blob_sha1=args.expected_git_blob_sha1,
        source_url=args.source_url,
        output_parquet=Path(args.output_parquet),
        output_json=Path(args.output_json),
    )
    if payload["status"] != "PASS":
        failed = [check["name"] for check in payload["checks"] if check["status"] == "FAIL"]
        print("sample metadata validation failed: " + ", ".join(failed))
        return 2
    print(
        "sample metadata validation passed: "
        f"{payload['join']['matched_sample_count']}/{payload['join']['plink_sample_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
