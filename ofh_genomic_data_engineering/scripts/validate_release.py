from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "variants": {"CHROM", "POS", "ID", "REF", "ALT"},
    "samples": {"IID"},
    "allele_frequencies": {"CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"},
    "variant_missingness": {"CHROM", "ID", "MISSING_CT", "OBS_CT", "F_MISS"},
    "sample_missingness": {"IID", "MISSING_CT", "OBS_CT", "F_MISS"},
    "hardy_weinberg": {"CHROM", "ID", "P"},
    "pca_scores": {"IID", "PC1"},
}

VARIANT_TABLES = {"allele_frequencies", "variant_missingness", "hardy_weinberg"}
SAMPLE_TABLES = {"sample_missingness", "pca_scores"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _canonical_sha256(payload: object) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _bounded_numeric(series: pd.Series, lower: float, upper: float) -> bool:
    values = pd.to_numeric(series, errors="coerce")
    return bool(values.notna().all() and values.between(lower, upper, inclusive="both").all())


def validate_release(
    *,
    source_inventory_path: Path,
    summary_path: Path,
    provenance_path: Path,
    parquet_dir: Path,
    bgen_path: Path,
    sample_path: Path,
    output_path: Path,
) -> dict[str, object]:
    inventory = json.loads(source_inventory_path.read_text())
    summary = json.loads(summary_path.read_text())
    provenance = json.loads(provenance_path.read_text())

    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, detail: object) -> None:
        checks.append({"name": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    parquet_files = summary.get("parquet_files", {})
    expected_tables = set(REQUIRED_COLUMNS)
    record(
        "table_inventory",
        set(parquet_files) == expected_tables,
        {"expected": sorted(expected_tables), "observed": sorted(parquet_files)},
    )

    tables: dict[str, pd.DataFrame] = {}
    for name in sorted(expected_tables):
        filename = parquet_files.get(name)
        path = parquet_dir / filename if filename else None
        exists = bool(path and path.is_file() and path.stat().st_size > 0)
        record(f"{name}.file_nonempty", exists, str(path) if path else "missing from summary")
        if not exists or path is None:
            continue
        frame = pd.read_parquet(path)
        tables[name] = frame
        missing = sorted(REQUIRED_COLUMNS[name] - set(frame.columns))
        record(f"{name}.required_columns", not missing, {"missing": missing})
        expected_rows = summary.get("row_counts", {}).get(name)
        record(
            f"{name}.row_count",
            expected_rows == len(frame),
            {"expected": expected_rows, "observed": len(frame)},
        )

    sample_count = summary.get("sample_count")
    variant_count = summary.get("variant_count")
    record(
        "sample_count_positive",
        isinstance(sample_count, int) and sample_count > 0,
        sample_count,
    )
    record(
        "variant_count_positive",
        isinstance(variant_count, int) and variant_count > 0,
        variant_count,
    )
    record(
        "source_sample_preservation",
        sample_count == inventory.get("sample_count"),
        {"source": inventory.get("sample_count"), "release": sample_count},
    )
    source_variant_count = inventory.get("variant_count")
    record(
        "variant_count_not_inflated",
        isinstance(source_variant_count, int)
        and isinstance(variant_count, int)
        and 0 < variant_count <= source_variant_count,
        {"source": source_variant_count, "release": variant_count},
    )

    if "samples" in tables and "IID" in tables["samples"]:
        sample_ids = tables["samples"]["IID"].astype(str)
        record("sample_ids_unique", sample_ids.is_unique, {"rows": len(sample_ids), "unique": sample_ids.nunique()})
        record(
            "samples_match_summary",
            len(sample_ids) == sample_count,
            {"rows": len(sample_ids), "summary": sample_count},
        )
        sample_set = set(sample_ids)
        for name in sorted(SAMPLE_TABLES):
            if name in tables and "IID" in tables[name]:
                other_ids = tables[name]["IID"].astype(str)
                record(
                    f"{name}.sample_ids_match",
                    set(other_ids) == sample_set and other_ids.is_unique,
                    {"rows": len(other_ids), "unique": other_ids.nunique()},
                )

    if "variants" in tables:
        variants = tables["variants"]
        key_columns = ["CHROM", "POS", "REF", "ALT"]
        if all(column in variants for column in key_columns):
            duplicated = variants.duplicated(key_columns).sum()
            record(
                "variant_keys_unique",
                duplicated == 0,
                {"rows": len(variants), "duplicate_keys": int(duplicated)},
            )
            record(
                "variants_match_summary",
                len(variants) == variant_count,
                {"rows": len(variants), "summary": variant_count},
            )

    for name in sorted(VARIANT_TABLES):
        if name in tables:
            record(
                f"{name}.variant_rows_match",
                len(tables[name]) == variant_count,
                {"rows": len(tables[name]), "summary": variant_count},
            )

    numeric_rules = [
        ("allele_frequencies", "ALT_FREQS"),
        ("variant_missingness", "F_MISS"),
        ("sample_missingness", "F_MISS"),
        ("hardy_weinberg", "P"),
    ]
    for table_name, column in numeric_rules:
        if table_name in tables and column in tables[table_name]:
            record(
                f"{table_name}.{column}_bounded_0_1",
                _bounded_numeric(tables[table_name][column], 0.0, 1.0),
                {"rows": len(tables[table_name])},
            )

    product_paths = [bgen_path, sample_path]
    product_paths.extend(
        parquet_dir / parquet_files[name]
        for name in sorted(expected_tables)
        if name in parquet_files
    )
    provenance_products = provenance.get("products", {})
    for path in product_paths:
        expected = provenance_products.get(path.name, {}).get("sha256")
        observed = sha256(path) if path.is_file() else None
        record(
            f"hash.{path.name}",
            isinstance(expected, str) and len(expected) == 64 and observed == expected,
            {"expected": expected, "observed": observed},
        )

    release_basis = {
        "source": provenance.get("source"),
        "normalised_vcf_sha256": inventory.get("normalised_vcf_sha256"),
        "parameters": provenance.get("parameters"),
        "summary": summary,
        "products": {
            path.name: provenance_products.get(path.name, {}).get("sha256")
            for path in product_paths
        },
    }
    release_id = _canonical_sha256(release_basis)
    passed = all(check["status"] == "PASS" for check in checks)
    payload = {
        "status": "PASS" if passed else "FAIL",
        "release_id": release_id,
        "source_inventory": inventory,
        "summary": summary,
        "checks": checks,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--source-inventory", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--provenance", required=True)
    p.add_argument("--parquet-dir", required=True)
    p.add_argument("--bgen", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--output", required=True)
    return p


def main() -> int:
    args = parser().parse_args()
    payload = validate_release(
        source_inventory_path=Path(args.source_inventory),
        summary_path=Path(args.summary),
        provenance_path=Path(args.provenance),
        parquet_dir=Path(args.parquet_dir),
        bgen_path=Path(args.bgen),
        sample_path=Path(args.sample),
        output_path=Path(args.output),
    )
    if payload["status"] != "PASS":
        failed = [check["name"] for check in payload["checks"] if check["status"] == "FAIL"]
        print("release validation failed: " + ", ".join(failed))
        return 2
    print(f"release validation passed: {payload['release_id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
