from pathlib import Path

import pandas as pd
import pytest

from scripts.validate_sample_metadata import git_blob_sha1, validate_sample_metadata


def _write_fixture(tmp_path: Path) -> dict[str, object]:
    pytest.importorskip("pyarrow")
    psam = tmp_path / "qc.psam"
    psam.write_text("#FID\tIID\tSEX\nS1_S1\tS1_S1\tNA\nS2_S2\tS2_S2\tNA\n")
    metadata = tmp_path / "panel.tsv"
    metadata.write_text(
        "sample\tpop\tsuper_pop\tgender\n"
        "S1\tGBR\tEUR\tmale\n"
        "S2\tFIN\tEUR\tfemale\n"
        "S3\tGBR\tEUR\tfemale\n"
    )
    return {
        "psam_path": psam,
        "metadata_path": metadata,
        "expected_git_blob_sha1": git_blob_sha1(metadata),
        "source_url": "https://example.org/panel.tsv",
        "output_parquet": tmp_path / "sample_metadata.parquet",
        "output_json": tmp_path / "metadata_validation.json",
    }


def test_sample_metadata_contract_passes_full_coverage(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    payload = validate_sample_metadata(**paths)
    assert payload["status"] == "PASS"
    assert payload["join"]["matched_sample_count"] == 2
    assert payload["join"]["coverage"] == 1.0
    assert len(payload["join"]["canonical_sample_ids_sha256"]) == 64
    assert len(payload["join"]["output_semantic_sha256"]) == 64
    frame = pd.read_parquet(paths["output_parquet"])
    assert frame["IID"].tolist() == ["S1_S1", "S2_S2"]
    assert frame["sample"].tolist() == ["S1", "S2"]
    assert frame["pop"].tolist() == ["GBR", "FIN"]


def test_sample_metadata_contract_rejects_missing_sample(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    paths["metadata_path"].write_text(
        "sample\tpop\tsuper_pop\tgender\nS1\tGBR\tEUR\tmale\n"
    )
    paths["expected_git_blob_sha1"] = git_blob_sha1(paths["metadata_path"])
    payload = validate_sample_metadata(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "metadata_full_coverage" in failed


def test_sample_metadata_contract_rejects_duplicate_metadata_id(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    paths["metadata_path"].write_text(
        "sample\tpop\tsuper_pop\tgender\n"
        "S1\tGBR\tEUR\tmale\n"
        "S1\tFIN\tEUR\tfemale\n"
        "S2\tFIN\tEUR\tfemale\n"
    )
    paths["expected_git_blob_sha1"] = git_blob_sha1(paths["metadata_path"])
    payload = validate_sample_metadata(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "metadata_sample_ids_unique" in failed


def test_sample_metadata_contract_rejects_noncanonical_plink_id(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    paths["psam_path"].write_text("#FID\tIID\tSEX\nS1_S1\tS1_OTHER\tNA\n")
    payload = validate_sample_metadata(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "plink_iid_canonical_form" in failed


def test_sample_metadata_contract_rejects_wrong_pinned_blob(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    paths["expected_git_blob_sha1"] = "0" * 40
    payload = validate_sample_metadata(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "metadata_git_blob_sha1" in failed
