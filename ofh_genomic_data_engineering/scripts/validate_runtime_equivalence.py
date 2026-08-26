from __future__ import annotations

import argparse
import json
from pathlib import Path


def _load(results: Path, relative: str) -> dict[str, object]:
    path = results / relative
    if not path.is_file():
        raise FileNotFoundError(f"missing runtime-equivalence input: {path}")
    return json.loads(path.read_text())


def _bgen_semantics(payload: dict[str, object]) -> dict[str, object]:
    keys = (
        "status",
        "contract",
        "allele_convention",
        "probability_bits",
        "frequency_tolerance",
        "sample_count",
        "variant_count",
        "sample_ids_sha256",
        "variant_identity_sha256",
    )
    return {key: payload.get(key) for key in keys}


def validate_runtime_equivalence(
    *,
    host_results: Path,
    candidate_results: Path,
    output_path: Path,
    expected_samples: int | None = None,
    expected_variants: int | None = None,
) -> dict[str, object]:
    host_delivery = _load(host_results, "00_source/delivery_validation.json")
    candidate_delivery = _load(candidate_results, "00_source/delivery_validation.json")
    host_bgen = _load(host_results, "05_bgen/bgen_validation.json")
    candidate_bgen = _load(candidate_results, "05_bgen/bgen_validation.json")
    host_summary = _load(host_results, "06_parquet/summary.json")
    candidate_summary = _load(candidate_results, "06_parquet/summary.json")
    host_query = _load(host_results, "06_parquet/query_validation.json")
    candidate_query = _load(candidate_results, "06_parquet/query_validation.json")
    host_release = _load(host_results, "08_release/release_validation.json")
    candidate_release = _load(candidate_results, "08_release/release_validation.json")

    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, detail: object) -> None:
        checks.append({"name": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    record(
        "delivery_status",
        host_delivery.get("status") == candidate_delivery.get("status") == "PASS",
        {"host": host_delivery.get("status"), "candidate": candidate_delivery.get("status")},
    )
    record(
        "delivery_action",
        host_delivery.get("action") == candidate_delivery.get("action") == "PROCESS",
        {"host": host_delivery.get("action"), "candidate": candidate_delivery.get("action")},
    )
    host_fingerprint = host_delivery.get("delivery", {}).get("delivery_fingerprint")
    candidate_fingerprint = candidate_delivery.get("delivery", {}).get("delivery_fingerprint")
    record(
        "delivery_fingerprint",
        isinstance(host_fingerprint, str)
        and len(host_fingerprint) == 64
        and host_fingerprint == candidate_fingerprint,
        {"host": host_fingerprint, "candidate": candidate_fingerprint},
    )

    host_samples = host_summary.get("sample_count")
    candidate_samples = candidate_summary.get("sample_count")
    record(
        "sample_count",
        host_samples == candidate_samples
        and (expected_samples is None or host_samples == expected_samples),
        {"host": host_samples, "candidate": candidate_samples, "expected": expected_samples},
    )
    host_variants = host_summary.get("variant_count")
    candidate_variants = candidate_summary.get("variant_count")
    record(
        "variant_count",
        host_variants == candidate_variants
        and (expected_variants is None or host_variants == expected_variants),
        {"host": host_variants, "candidate": candidate_variants, "expected": expected_variants},
    )

    host_semantic = host_summary.get("semantic_hashes")
    candidate_semantic = candidate_summary.get("semantic_hashes")
    record(
        "semantic_hashes",
        isinstance(host_semantic, dict)
        and len(host_semantic) == 7
        and host_semantic == candidate_semantic,
        {"host": host_semantic, "candidate": candidate_semantic},
    )

    host_bgen_semantics = _bgen_semantics(host_bgen)
    candidate_bgen_semantics = _bgen_semantics(candidate_bgen)
    record(
        "bgen_roundtrip_semantics",
        host_bgen_semantics == candidate_bgen_semantics
        and host_bgen.get("status") == candidate_bgen.get("status") == "PASS",
        {"host": host_bgen_semantics, "candidate": candidate_bgen_semantics},
    )
    record(
        "bgen_roundtrip_frequency_bound",
        isinstance(host_bgen.get("max_abs_alt_frequency_diff"), (int, float))
        and isinstance(candidate_bgen.get("max_abs_alt_frequency_diff"), (int, float))
        and host_bgen.get("max_abs_alt_frequency_diff") <= host_bgen.get("frequency_tolerance")
        and candidate_bgen.get("max_abs_alt_frequency_diff") <= candidate_bgen.get("frequency_tolerance"),
        {
            "host_max_abs_diff": host_bgen.get("max_abs_alt_frequency_diff"),
            "candidate_max_abs_diff": candidate_bgen.get("max_abs_alt_frequency_diff"),
            "tolerance": host_bgen.get("frequency_tolerance"),
        },
    )

    record(
        "query_status",
        host_query.get("status") == candidate_query.get("status") == "PASS",
        {"host": host_query.get("status"), "candidate": candidate_query.get("status")},
    )
    record(
        "query_result",
        host_query.get("total_variants") == candidate_query.get("total_variants")
        and host_query.get("region_query") == candidate_query.get("region_query"),
        {
            "host_total_variants": host_query.get("total_variants"),
            "candidate_total_variants": candidate_query.get("total_variants"),
            "host_region": host_query.get("region_query"),
            "candidate_region": candidate_query.get("region_query"),
        },
    )

    host_identity_version = host_release.get("release_identity", {}).get("version")
    candidate_identity_version = candidate_release.get("release_identity", {}).get("version")
    record(
        "release_identity_version",
        host_identity_version == candidate_identity_version == 3,
        {"host": host_identity_version, "candidate": candidate_identity_version},
    )
    record(
        "release_status",
        host_release.get("status") == candidate_release.get("status") == "PASS",
        {"host": host_release.get("status"), "candidate": candidate_release.get("status")},
    )
    host_release_id = host_release.get("release_id")
    candidate_release_id = candidate_release.get("release_id")
    record(
        "semantic_release_id",
        isinstance(host_release_id, str)
        and len(host_release_id) == 64
        and host_release_id == candidate_release_id,
        {"host": host_release_id, "candidate": candidate_release_id},
    )

    passed = all(check["status"] == "PASS" for check in checks)
    payload = {
        "status": "PASS" if passed else "FAIL",
        "contract": "cross-runtime-semantic-equivalence",
        "release_identity_version": host_identity_version,
        "release_id": host_release_id if host_release_id == candidate_release_id else None,
        "delivery_fingerprint": host_fingerprint if host_fingerprint == candidate_fingerprint else None,
        "sample_count": host_samples if host_samples == candidate_samples else None,
        "variant_count": host_variants if host_variants == candidate_variants else None,
        "semantic_hashes": host_semantic if host_semantic == candidate_semantic else None,
        "bgen_roundtrip": host_bgen_semantics if host_bgen_semantics == candidate_bgen_semantics else None,
        "checks": checks,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--host-results", required=True)
    p.add_argument("--candidate-results", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--expected-samples", type=int)
    p.add_argument("--expected-variants", type=int)
    return p


def main() -> int:
    args = parser().parse_args()
    payload = validate_runtime_equivalence(
        host_results=Path(args.host_results),
        candidate_results=Path(args.candidate_results),
        output_path=Path(args.output),
        expected_samples=args.expected_samples,
        expected_variants=args.expected_variants,
    )
    if payload["status"] != "PASS":
        failed = [check["name"] for check in payload["checks"] if check["status"] == "FAIL"]
        print("runtime equivalence failed: " + ", ".join(failed))
        return 2
    print(f"runtime equivalence passed: {payload['release_id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
