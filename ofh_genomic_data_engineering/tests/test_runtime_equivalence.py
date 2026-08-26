import json
from pathlib import Path

from scripts.validate_runtime_equivalence import validate_runtime_equivalence


def _write_runtime(results: Path, *, semantic_suffix: str = "a") -> None:
    (results / "00_source").mkdir(parents=True)
    (results / "06_parquet").mkdir(parents=True)
    (results / "08_release").mkdir(parents=True)
    fingerprint = "f" * 64
    semantic_hashes = {name: semantic_suffix * 64 for name in ("variants", "samples", "allele_frequencies", "variant_missingness", "sample_missingness", "hardy_weinberg", "pca_scores")}
    (results / "00_source/delivery_validation.json").write_text(json.dumps({"status": "PASS", "action": "PROCESS", "delivery": {"delivery_fingerprint": fingerprint}}))
    (results / "06_parquet/summary.json").write_text(json.dumps({"sample_count": 90, "variant_count": 127171, "semantic_hashes": semantic_hashes}))
    (results / "06_parquet/query_validation.json").write_text(json.dumps({"status": "PASS", "total_variants": 127171, "region_query": {"start": 1, "end": 2, "matched_variants": 10}}))
    release_id = "e" * 64 if semantic_suffix == "a" else "d" * 64
    (results / "08_release/release_validation.json").write_text(json.dumps({"status": "PASS", "release_id": release_id, "release_identity": {"version": 2}}))


def test_runtime_equivalence_passes_for_matching_semantics(tmp_path: Path) -> None:
    host = tmp_path / "host"
    candidate = tmp_path / "candidate"
    _write_runtime(host)
    _write_runtime(candidate)
    payload = validate_runtime_equivalence(host_results=host, candidate_results=candidate, output_path=tmp_path / "evidence.json", expected_samples=90, expected_variants=127171)
    assert payload["status"] == "PASS"
    assert payload["release_id"] == "e" * 64
    assert all(check["status"] == "PASS" for check in payload["checks"])


def test_runtime_equivalence_rejects_semantic_drift(tmp_path: Path) -> None:
    host = tmp_path / "host"
    candidate = tmp_path / "candidate"
    _write_runtime(host)
    _write_runtime(candidate, semantic_suffix="b")
    payload = validate_runtime_equivalence(host_results=host, candidate_results=candidate, output_path=tmp_path / "evidence.json")
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "semantic_hashes" in failed
    assert "semantic_release_id" in failed
