from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def _read_plink_table(path: Path) -> pd.DataFrame:
    skiprows = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                skiprows += 1
                continue
            break
    frame = pd.read_csv(path, sep=r"\s+", skiprows=skiprows)
    frame.columns = [str(column).lstrip("#") for column in frame.columns]
    return frame


def _sha256_lines(lines: list[str]) -> str:
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _normalise_variant_table(frame: pd.DataFrame) -> pd.DataFrame:
    required = ["CHROM", "POS", "ID", "REF", "ALT"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"missing PVAR columns: {missing}")
    out = frame[required].copy()
    out["CHROM"] = out["CHROM"].astype(str)
    out["POS"] = pd.to_numeric(out["POS"], errors="raise").astype("int64")
    for column in ("ID", "REF", "ALT"):
        out[column] = out[column].astype(str)
    return out


def _variant_hash(frame: pd.DataFrame) -> str:
    return _sha256_lines(
        [
            "\t".join((str(row.CHROM), str(row.POS), row.ID, row.REF, row.ALT))
            for row in frame.itertuples(index=False)
        ]
    )


def _sample_hash(frame: pd.DataFrame) -> str:
    if "IID" not in frame.columns:
        raise ValueError("PSAM is missing IID")
    return _sha256_lines(frame["IID"].astype(str).tolist())


def validate_bgen_roundtrip(
    *,
    source_pvar: Path,
    source_psam: Path,
    source_afreq: Path,
    roundtrip_pvar: Path,
    roundtrip_psam: Path,
    roundtrip_afreq: Path,
    output_path: Path,
    frequency_tolerance: float,
    probability_bits: int = 16,
) -> dict[str, object]:
    if frequency_tolerance < 0:
        raise ValueError("frequency tolerance must be non-negative")

    source_variants = _normalise_variant_table(_read_plink_table(source_pvar))
    roundtrip_variants = _normalise_variant_table(_read_plink_table(roundtrip_pvar))
    source_samples = _read_plink_table(source_psam)
    roundtrip_samples = _read_plink_table(roundtrip_psam)
    source_freq = _read_plink_table(source_afreq)
    roundtrip_freq = _read_plink_table(roundtrip_afreq)

    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, detail: object) -> None:
        checks.append({"name": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    source_iids = source_samples["IID"].astype(str).tolist() if "IID" in source_samples else []
    roundtrip_iids = roundtrip_samples["IID"].astype(str).tolist() if "IID" in roundtrip_samples else []
    source_sample_hash = _sample_hash(source_samples)
    roundtrip_sample_hash = _sample_hash(roundtrip_samples)
    record(
        "sample_order_identity",
        source_iids == roundtrip_iids,
        {"source_count": len(source_iids), "roundtrip_count": len(roundtrip_iids)},
    )
    record(
        "sample_identity_hash",
        source_sample_hash == roundtrip_sample_hash,
        {"source": source_sample_hash, "roundtrip": roundtrip_sample_hash},
    )

    source_variant_hash = _variant_hash(source_variants)
    roundtrip_variant_hash = _variant_hash(roundtrip_variants)
    record(
        "variant_count_identity",
        len(source_variants) == len(roundtrip_variants),
        {"source": len(source_variants), "roundtrip": len(roundtrip_variants)},
    )
    record(
        "variant_order_allele_identity",
        source_variants.equals(roundtrip_variants),
        {"source_hash": source_variant_hash, "roundtrip_hash": roundtrip_variant_hash},
    )

    required_freq = ["CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"]
    missing_source = [column for column in required_freq if column not in source_freq.columns]
    missing_roundtrip = [column for column in required_freq if column not in roundtrip_freq.columns]
    record(
        "frequency_columns",
        not missing_source and not missing_roundtrip,
        {"source_missing": missing_source, "roundtrip_missing": missing_roundtrip},
    )

    max_abs_diff: float | None = None
    mean_abs_diff: float | None = None
    obs_ct_mismatches: int | None = None
    frequency_order_identity = False
    frequency_within_tolerance = False

    if not missing_source and not missing_roundtrip and len(source_freq) == len(roundtrip_freq):
        source_keys = source_freq[["CHROM", "ID", "REF", "ALT"]].astype(str).reset_index(drop=True)
        roundtrip_keys = roundtrip_freq[["CHROM", "ID", "REF", "ALT"]].astype(str).reset_index(drop=True)
        frequency_order_identity = source_keys.equals(roundtrip_keys)

        source_obs = pd.to_numeric(source_freq["OBS_CT"], errors="coerce")
        roundtrip_obs = pd.to_numeric(roundtrip_freq["OBS_CT"], errors="coerce")
        obs_ct_mismatches = int((source_obs != roundtrip_obs).sum())

        source_alt = pd.to_numeric(source_freq["ALT_FREQS"], errors="coerce")
        roundtrip_alt = pd.to_numeric(roundtrip_freq["ALT_FREQS"], errors="coerce")
        finite = source_alt.notna() & roundtrip_alt.notna()
        if bool(finite.all()) and len(source_alt) > 0:
            diffs = (source_alt - roundtrip_alt).abs()
            max_abs_diff = float(diffs.max())
            mean_abs_diff = float(diffs.mean())
            frequency_within_tolerance = bool((diffs <= frequency_tolerance).all())

    record(
        "frequency_row_count_identity",
        len(source_freq) == len(roundtrip_freq),
        {"source": len(source_freq), "roundtrip": len(roundtrip_freq)},
    )
    record(
        "frequency_variant_order_identity",
        frequency_order_identity,
        {"rows": len(source_freq)},
    )
    record(
        "frequency_observation_count_identity",
        obs_ct_mismatches == 0,
        {"mismatches": obs_ct_mismatches},
    )
    record(
        "alt_frequency_tolerance",
        frequency_within_tolerance,
        {
            "tolerance": frequency_tolerance,
            "max_abs_diff": max_abs_diff,
            "mean_abs_diff": mean_abs_diff,
        },
    )

    passed = all(check["status"] == "PASS" for check in checks)
    payload = {
        "status": "PASS" if passed else "FAIL",
        "contract": "bgen-1.2-roundtrip",
        "allele_convention": "ref-first",
        "probability_bits": probability_bits,
        "frequency_tolerance": frequency_tolerance,
        "sample_count": len(source_iids),
        "variant_count": len(source_variants),
        "sample_ids_sha256": source_sample_hash if source_sample_hash == roundtrip_sample_hash else None,
        "variant_identity_sha256": source_variant_hash if source_variant_hash == roundtrip_variant_hash else None,
        "max_abs_alt_frequency_diff": max_abs_diff,
        "mean_abs_alt_frequency_diff": mean_abs_diff,
        "obs_ct_mismatch_count": obs_ct_mismatches,
        "checks": checks,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--source-pvar", required=True)
    p.add_argument("--source-psam", required=True)
    p.add_argument("--source-afreq", required=True)
    p.add_argument("--roundtrip-pvar", required=True)
    p.add_argument("--roundtrip-psam", required=True)
    p.add_argument("--roundtrip-afreq", required=True)
    p.add_argument("--frequency-tolerance", type=float, default=1e-4)
    p.add_argument("--probability-bits", type=int, default=16)
    p.add_argument("--output", required=True)
    return p


def main() -> int:
    args = parser().parse_args()
    payload = validate_bgen_roundtrip(
        source_pvar=Path(args.source_pvar),
        source_psam=Path(args.source_psam),
        source_afreq=Path(args.source_afreq),
        roundtrip_pvar=Path(args.roundtrip_pvar),
        roundtrip_psam=Path(args.roundtrip_psam),
        roundtrip_afreq=Path(args.roundtrip_afreq),
        output_path=Path(args.output),
        frequency_tolerance=args.frequency_tolerance,
        probability_bits=args.probability_bits,
    )
    if payload["status"] != "PASS":
        failed = [check["name"] for check in payload["checks"] if check["status"] == "FAIL"]
        print("BGEN round-trip validation failed: " + ", ".join(failed))
        return 2
    print(
        "BGEN round-trip validation passed: "
        f"{payload['sample_count']} samples, {payload['variant_count']} variants, "
        f"max ALT-frequency diff={payload['max_abs_alt_frequency_diff']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
