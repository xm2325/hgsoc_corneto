import json
from pathlib import Path

from scripts.validate_runtime_equivalence import validate_runtime_equivalence


def _write_runtime(
    results: Path,
    *,
    semantic_suffix: str = "a",
    bgen_suffix: str = "c",
    metadata_suffix: str = "1",
    plink_seed: int = 20260826,
) -> None:
    (results / "00_source").mkdir(parents=True)
    (results / "05_bgen").mkdir(parents=True)
    (results / "05_metadata").mkdir(parents=True)
    (results / "06_parquet").mkdir(parents=True)
    (results / "08_release").mkdir(parents=True)
    fingerprint = "f" * 64
    semantic_hashes = {
        name: semantic_suffix * 64
        for name in (
            "variants",
            "samples",
            "allele_frequencies",
            "variant_missingness",
            "sample_missingness",
            "hardy_weinberg",
            "pca_scores",
        )
    }
    (results / "00_source/delivery_validation.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "action": "PROCESS",
                "delivery": {"delivery_fingerprint": fingerprint},
            }
        )
    )
    (results / "05_bgen/bgen_validation.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "contract": "bgen-1.2-roundtrip",
                "allele_convention": "ref-first",
                "probability_bits": 16,
                "frequency_tolerance": 0.0001,
                "sample_count": 90,
                "variant_count": 127171,
                "sample_ids_sha256": bgen_suffix * 64,
                "variant_identity_sha256": "d" * 64,
                "max_abs_alt_frequency_diff": 0.0,
            }
        )
    )
    metadata = {
        "status": "PASS",
        "contract": "sample-metadata-join-v1",
        "source": {
            "git_blob_sha1": "b" * 40,
            "sha256": "2" * 64,
            "row_count": 2504,
        },
        "join": {
            "plink_sample_count": 90,
            "matched_sample_count": 90,
            "coverage": 1.0,
            "canonical_sample_ids_sha256": "3" * 64,
            "output_semantic_sha256": metadata_suffix * 64,
        },
    }
    (results / "05_metadata/metadata_validation.json").write_text(json.dumps(metadata))
    (results / "06_parquet/summary.json").write_text(
        json.dumps(
            {
                "sample_count": 90,
                "variant_count": 127171,
                "semantic_hashes": semantic_hashes,
            }
        )
    )
    (results / "06_parquet/query_validation.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "total_variants": 127171,
                "region_query": {"start": 1, "end": 2, "matched_variants": 10},
            }
        )
    )
    release_id = (
        "e" * 64
        if semantic_suffix == "a"
        and bgen_suffix == "c"
        and metadata_suffix == "1"
        and plink_seed == 20260826
        else "d" * 64
    )
    (results / "08_release/release_validation.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "release_id": release_id,
                "release_identity": {
                    "version": 4,
                    "basis": {
                        "parameters": {
                            "plink_seed": plink_seed,
                            "plink_threads": 2,
                            "plink_memory_mb": 3000,
                        }
                    },
                },
            }
        )
    )


def test_runtime_equivalence_passes_for_matching_semantics(tmp_path: Path) -> None:
    host = tmp_path / "host"
    candidate = tmp_path / "candidate"
    _write_runtime(host)
    _write_runtime(candidate)
    payload = validate_runtime_equivalence(
        host_results=host,
        candidate_results=candidate,
        output_path=tmp_path / "evidence.json",
        expected_samples=90,
        expected_variants=127171,
    )
    assert payload["status"] == "PASS"
    assert payload["release_id"] == "e" * 64
    assert payload["release_identity_version"] == 4
    assert payload["bgen_roundtrip"]["allele_convention"] == "ref-first"
    assert payload["sample_metadata_join"]["coverage"] == 1.0
    assert payload["pca_execution_parameters"] == {
        "plink_seed": 20260826,
        "plink_threads": 2,
        "plink_memory_mb": 3000,
    }
    assert all(check["status"] == "PASS" for check in payload["checks"])


def test_runtime_equivalence_rejects_parquet_semantic_drift(tmp_path: Path) -> None:
    host = tmp_path / "host"
    candidate = tmp_path / "candidate"
    _write_runtime(host)
    _write_runtime(candidate, semantic_suffix="b")
    payload = validate_runtime_equivalence(
        host_results=host,
        candidate_results=candidate,
        output_path=tmp_path / "evidence.json",
    )
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "semantic_hashes" in failed
    assert "semantic_release_id" in failed


def test_runtime_equivalence_rejects_bgen_semantic_drift(tmp_path: Path) -> None:
    host = tmp_path / "host"
    candidate = tmp_path / "candidate"
    _write_runtime(host)
    _write_runtime(candidate, bgen_suffix="b")
    payload = validate_runtime_equivalence(
        host_results=host,
        candidate_results=candidate,
        output_path=tmp_path / "evidence.json",
    )
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "bgen_roundtrip_semantics" in failed
    assert "semantic_release_id" in failed


def test_runtime_equivalence_rejects_metadata_semantic_drift(tmp_path: Path) -> None:
    host = tmp_path / "host"
    candidate = tmp_path / "candidate"
    _write_runtime(host)
    _write_runtime(candidate, metadata_suffix="4")
    payload = validate_runtime_equivalence(
        host_results=host,
        candidate_results=candidate,
        output_path=tmp_path / "evidence.json",
    )
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "sample_metadata_semantics" in failed
    assert "semantic_release_id" in failed


def test_runtime_equivalence_rejects_pca_execution_parameter_drift(tmp_path: Path) -> None:
    host = tmp_path / "host"
    candidate = tmp_path / "candidate"
    _write_runtime(host)
    _write_runtime(candidate, plink_seed=20260827)
    payload = validate_runtime_equivalence(
        host_results=host,
        candidate_results=candidate,
        output_path=tmp_path / "evidence.json",
    )
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "pca_execution_parameters" in failed
    assert "semantic_release_id" in failed
