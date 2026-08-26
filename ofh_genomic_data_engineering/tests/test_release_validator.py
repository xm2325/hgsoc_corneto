import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.metrics_to_parquet import semantic_table_sha256
from scripts.validate_release import sha256, validate_release


def _write_release(tmp_path: Path) -> dict[str, Path]:
    pytest.importorskip("pyarrow")
    parquet_dir = tmp_path / "parquet"
    parquet_dir.mkdir()

    tables = {
        "variants": pd.DataFrame({"CHROM": ["22", "22"], "POS": ["1", "2"], "ID": ["rs1", "rs2"], "REF": ["A", "C"], "ALT": ["G", "T"]}),
        "samples": pd.DataFrame({"FID": ["S1", "S2"], "IID": ["S1", "S2"]}),
        "allele_frequencies": pd.DataFrame({"CHROM": ["22", "22"], "ID": ["rs1", "rs2"], "REF": ["A", "C"], "ALT": ["G", "T"], "ALT_FREQS": ["0.25", "0.5"], "OBS_CT": ["4", "4"]}),
        "variant_missingness": pd.DataFrame({"CHROM": ["22", "22"], "ID": ["rs1", "rs2"], "MISSING_CT": ["0", "0"], "OBS_CT": ["2", "2"], "F_MISS": ["0", "0"]}),
        "sample_missingness": pd.DataFrame({"FID": ["S1", "S2"], "IID": ["S1", "S2"], "MISSING_CT": ["0", "0"], "OBS_CT": ["2", "2"], "F_MISS": ["0", "0"]}),
        "hardy_weinberg": pd.DataFrame({"CHROM": ["22", "22"], "ID": ["rs1", "rs2"], "P": ["1", "0.75"]}),
        "pca_scores": pd.DataFrame({"FID": ["S1", "S2"], "IID": ["S1", "S2"], "PC1": ["0.1", "-0.1"]}),
    }

    parquet_files = {}
    row_counts = {}
    semantic_hashes = {}
    for name, frame in tables.items():
        path = parquet_dir / f"{name}.parquet"
        frame.to_parquet(path, index=False)
        stored = pd.read_parquet(path)
        parquet_files[name] = path.name
        row_counts[name] = len(stored)
        semantic_hashes[name] = semantic_table_sha256(stored)

    summary = {"sample_count": 2, "variant_count": 2, "row_counts": row_counts, "parquet_files": parquet_files, "semantic_hashes": semantic_hashes}
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(summary))

    inventory = {"sample_count": 2, "sample_ids_sha256": "a" * 64, "variant_count": 3, "normalised_vcf_sha256": "b" * 64}
    inventory_path = tmp_path / "source_inventory.json"
    inventory_path.write_text(json.dumps(inventory))

    bgen_path = tmp_path / "analysis_ready.bgen"
    bgen_path.write_bytes(b"BGEN fixture")
    sample_path = tmp_path / "analysis_ready.sample"
    sample_path.write_text("ID_1 ID_2\n0 0\nS1 S1\nS2 S2\n")
    bgen_validation = {
        "status": "PASS",
        "contract": "bgen-1.2-roundtrip",
        "allele_convention": "ref-first",
        "probability_bits": 16,
        "frequency_tolerance": 0.0001,
        "sample_count": 2,
        "variant_count": 2,
        "sample_ids_sha256": "e" * 64,
        "variant_identity_sha256": "f" * 64,
        "max_abs_alt_frequency_diff": 0.0,
        "mean_abs_alt_frequency_diff": 0.0,
        "obs_ct_mismatch_count": 0,
        "checks": [{"name": "fixture", "status": "PASS", "detail": {}}],
    }
    bgen_validation_path = tmp_path / "bgen_validation.json"
    bgen_validation_path.write_text(json.dumps(bgen_validation))

    product_paths = [bgen_path, sample_path, bgen_validation_path, *sorted(parquet_dir.glob("*.parquet"))]
    provenance = {
        "source": {"url": "https://example.org/input.vcf.gz", "sha256": "c" * 64},
        "delivery": {"delivery_fingerprint": "d" * 64, "reference_genome": "GRCh37"},
        "parameters": {"geno": "0.02", "maf": "0.01", "hwe": "1e-6", "delivery_fingerprint": "d" * 64, "reference_genome": "GRCh37"},
        "bgen_roundtrip": bgen_validation,
        "products": {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in product_paths},
    }
    provenance_path = tmp_path / "provenance.json"
    provenance_path.write_text(json.dumps(provenance))

    return {
        "source_inventory_path": inventory_path,
        "summary_path": summary_path,
        "provenance_path": provenance_path,
        "parquet_dir": parquet_dir,
        "bgen_path": bgen_path,
        "sample_path": sample_path,
        "bgen_validation_path": bgen_validation_path,
        "output_path": tmp_path / "release_validation.json",
    }


def _refresh_product_hash(paths: dict[str, Path], changed: Path) -> None:
    provenance = json.loads(paths["provenance_path"].read_text())
    provenance["products"][changed.name] = {"bytes": changed.stat().st_size, "sha256": sha256(changed)}
    paths["provenance_path"].write_text(json.dumps(provenance))


def test_release_contract_passes_for_consistent_products(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    payload = validate_release(**paths)
    assert payload["status"] == "PASS"
    assert len(payload["release_id"]) == 64
    assert payload["release_identity"]["version"] == 3
    assert payload["release_identity"]["basis"]["bgen_contract"]["allele_convention"] == "ref-first"
    assert all(check["status"] == "PASS" for check in payload["checks"])


def test_release_contract_rejects_duplicate_sample_ids(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    samples = pd.read_parquet(paths["parquet_dir"] / "samples.parquet")
    samples["IID"] = ["S1", "S1"]
    changed = paths["parquet_dir"] / "samples.parquet"
    samples.to_parquet(changed, index=False)
    _refresh_product_hash(paths, changed)
    payload = validate_release(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "sample_ids_unique" in failed
    assert "samples.semantic_hash" in failed


def test_release_contract_rejects_hash_mismatch(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    provenance = json.loads(paths["provenance_path"].read_text())
    provenance["products"]["analysis_ready.bgen"]["sha256"] = "0" * 64
    paths["provenance_path"].write_text(json.dumps(provenance))
    payload = validate_release(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "hash.analysis_ready.bgen" in failed


def test_release_contract_rejects_declared_semantic_hash_mismatch(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    summary = json.loads(paths["summary_path"].read_text())
    summary["semantic_hashes"]["variants"] = "0" * 64
    paths["summary_path"].write_text(json.dumps(summary))
    payload = validate_release(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "variants.semantic_hash" in failed


def test_release_contract_rejects_failed_bgen_roundtrip(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    bgen = json.loads(paths["bgen_validation_path"].read_text())
    bgen["status"] = "FAIL"
    bgen["max_abs_alt_frequency_diff"] = 0.01
    paths["bgen_validation_path"].write_text(json.dumps(bgen))
    _refresh_product_hash(paths, paths["bgen_validation_path"])
    payload = validate_release(**paths)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "bgen_roundtrip.status" in failed
    assert "bgen_roundtrip.frequency_tolerance" in failed
    assert "bgen_roundtrip.provenance_binding" in failed


def test_release_id_is_independent_of_parquet_serialisation(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    first = validate_release(**paths)
    assert first["status"] == "PASS"
    changed = paths["parquet_dir"] / "variants.parquet"
    variants = pd.read_parquet(changed)
    before_file_sha = sha256(changed)
    variants.to_parquet(changed, index=False, compression="gzip")
    assert sha256(changed) != before_file_sha
    _refresh_product_hash(paths, changed)
    second = validate_release(**paths)
    assert second["status"] == "PASS"
    assert second["release_id"] == first["release_id"]


def test_release_id_changes_when_semantic_content_changes(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    first = validate_release(**paths)
    assert first["status"] == "PASS"
    changed = paths["parquet_dir"] / "pca_scores.parquet"
    pca = pd.read_parquet(changed)
    pca.loc[0, "PC1"] = "0.2"
    pca.to_parquet(changed, index=False)
    _refresh_product_hash(paths, changed)
    summary = json.loads(paths["summary_path"].read_text())
    summary["semantic_hashes"]["pca_scores"] = semantic_table_sha256(pd.read_parquet(changed))
    paths["summary_path"].write_text(json.dumps(summary))
    second = validate_release(**paths)
    assert second["status"] == "PASS"
    assert second["release_id"] != first["release_id"]


def test_release_id_changes_when_bgen_contract_changes(tmp_path: Path) -> None:
    paths = _write_release(tmp_path)
    first = validate_release(**paths)
    assert first["status"] == "PASS"
    bgen = json.loads(paths["bgen_validation_path"].read_text())
    bgen["frequency_tolerance"] = 0.0002
    paths["bgen_validation_path"].write_text(json.dumps(bgen))
    _refresh_product_hash(paths, paths["bgen_validation_path"])
    provenance = json.loads(paths["provenance_path"].read_text())
    provenance["bgen_roundtrip"] = bgen
    paths["provenance_path"].write_text(json.dumps(provenance))
    second = validate_release(**paths)
    assert second["status"] == "PASS"
    assert second["release_id"] != first["release_id"]
