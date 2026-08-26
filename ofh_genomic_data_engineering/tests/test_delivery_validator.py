from __future__ import annotations

import gzip
import json
from pathlib import Path

from scripts.validate_delivery import sample_ids_sha256, sha256, validate_delivery


def write_vcf(path: Path, sample_ids: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids) + "\n")
        fh.write("22\t16050075\trs1\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(["0/1"] * len(sample_ids)) + "\n")


def write_manifest(path: Path, source: Path, sample_ids: list[str], **overrides: object) -> dict[str, object]:
    manifest: dict[str, object] = {
        "schema_version": 1,
        "delivery_id": "provider-2026-08-26-batch-001",
        "provider": "example-provider",
        "source_uri": "https://example.org/batch.vcf.gz",
        "source_format": "vcf.gz",
        "reference_genome": "GRCh37",
        "sha256": sha256(source),
        "sample_roster": {"count": len(sample_ids), "ids_sha256": sample_ids_sha256(sample_ids)},
    }
    manifest.update(overrides)
    path.write_text(json.dumps(manifest))
    return manifest


def test_valid_delivery_is_accepted_for_processing(tmp_path: Path) -> None:
    source = tmp_path / "batch.vcf.gz"; samples = ["S02", "S01"]; write_vcf(source, samples)
    manifest = tmp_path / "manifest.json"; write_manifest(manifest, source, samples)
    result = validate_delivery(manifest_path=manifest, source_path=source, output_path=tmp_path / "validation.json", required_reference_genome="GRCh37", expected_source_uri="https://example.org/batch.vcf.gz")
    assert result["status"] == "PASS" and result["action"] == "PROCESS" and result["should_process"] is True
    assert len(result["delivery"]["delivery_fingerprint"]) == 64


def test_checksum_mismatch_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "batch.vcf.gz"; samples = ["S01"]; write_vcf(source, samples)
    manifest = tmp_path / "manifest.json"; write_manifest(manifest, source, samples, sha256="0" * 64)
    result = validate_delivery(manifest_path=manifest, source_path=source, output_path=tmp_path / "validation.json", required_reference_genome="GRCh37")
    assert result["status"] == "FAIL" and result["action"] == "REJECT"
    assert "source.sha256" in {c["name"] for c in result["checks"] if c["status"] == "FAIL"}


def test_wrong_reference_genome_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "batch.vcf.gz"; samples = ["S01"]; write_vcf(source, samples)
    manifest = tmp_path / "manifest.json"; write_manifest(manifest, source, samples, reference_genome="GRCh38")
    result = validate_delivery(manifest_path=manifest, source_path=source, output_path=tmp_path / "validation.json", required_reference_genome="GRCh37")
    assert result["status"] == "FAIL"
    assert "manifest.reference_genome" in {c["name"] for c in result["checks"] if c["status"] == "FAIL"}


def test_sample_roster_mismatch_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "batch.vcf.gz"; samples = ["S01", "S02"]; write_vcf(source, samples)
    manifest = tmp_path / "manifest.json"; write_manifest(manifest, source, samples, sample_roster={"count": 2, "ids_sha256": "f" * 64})
    result = validate_delivery(manifest_path=manifest, source_path=source, output_path=tmp_path / "validation.json", required_reference_genome="GRCh37")
    assert result["status"] == "FAIL"
    assert "source.sample_roster_sha256" in {c["name"] for c in result["checks"] if c["status"] == "FAIL"}


def test_exact_duplicate_delivery_becomes_noop(tmp_path: Path) -> None:
    source = tmp_path / "batch.vcf.gz"; samples = ["S01", "S02"]; write_vcf(source, samples)
    manifest = tmp_path / "manifest.json"; write_manifest(manifest, source, samples)
    first = validate_delivery(manifest_path=manifest, source_path=source, output_path=tmp_path / "first.json", required_reference_genome="GRCh37")
    registry = tmp_path / "registry.json"; registry.write_text(json.dumps({"deliveries": [{"delivery_id": first["delivery"]["delivery_id"], "source_sha256": first["source_observed"]["sha256"], "delivery_fingerprint": first["delivery"]["delivery_fingerprint"]}]}))
    duplicate = validate_delivery(manifest_path=manifest, source_path=source, output_path=tmp_path / "duplicate.json", required_reference_genome="GRCh37", registry_path=registry)
    assert duplicate["status"] == "PASS" and duplicate["action"] == "NOOP" and duplicate["should_process"] is False


def test_delivery_id_collision_with_changed_content_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "batch.vcf.gz"; samples = ["S01"]; write_vcf(source, samples)
    manifest = tmp_path / "manifest.json"; write_manifest(manifest, source, samples)
    registry = tmp_path / "registry.json"; registry.write_text(json.dumps({"deliveries": [{"delivery_id": "provider-2026-08-26-batch-001", "source_sha256": "0" * 64, "delivery_fingerprint": "1" * 64}]}))
    result = validate_delivery(manifest_path=manifest, source_path=source, output_path=tmp_path / "validation.json", required_reference_genome="GRCh37", registry_path=registry)
    assert result["status"] == "FAIL" and result["action"] == "REJECT" and result["should_process"] is False
    collision = next(c for c in result["checks"] if c["name"] == "registry.delivery_id_collision")
    assert collision["status"] == "FAIL"
