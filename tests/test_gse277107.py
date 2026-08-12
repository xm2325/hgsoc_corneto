import csv
import gzip
import importlib.util
import json
import math
import subprocess
from pathlib import Path

import pytest

from hgsoc_corneto.external.gse277107 import (
    GSE277107Error,
    audit_prepared_dataset,
    build_paired_metadata,
    prepare_dataset,
)
from hgsoc_corneto.io import sha256

ROOT = Path(__file__).resolve().parents[1]
REGULATORY_SCRIPT = ROOT / "scripts/run_corneto_regulatory_pilot.py"
SPEC = importlib.util.spec_from_file_location("gse277107_regulatory_loader", REGULATORY_SCRIPT)
REGULATORY_LOADER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(REGULATORY_LOADER)


def _write_soft(path: Path, *, omit_omentum: bool = False) -> None:
    records = [
        ("GSM1", "156A_OV_HGSC_ovary", "WM156A_OV", "ovary", "SRX1"),
        ("GSM2", "156C_OM_HGSC_omentum", "WM156C_OM", "omentum", "SRX2"),
        ("GSM3", "1859A_OV_HGSC_ovary", "WM1859A_OV", "ovary", "SRX3"),
        ("GSM4", "1859B_OM_HGSCS_omentum", "WM1859B_OM", "omentum", "SRX4"),
    ]
    if omit_omentum:
        records.pop()
    lines = [
        "^SERIES = GSE277107",
        "!Series_geo_accession = GSE277107",
        "!Series_overall_design = matched primary site (ovary / adexa) and a common "
        "secondary site (omentum)",
    ]
    for gsm, title, description, tissue, srx in records:
        source = f"High Grade Serous Ovarian Cancer_{tissue}"
        lines.extend(
            [
                f"^SAMPLE = {gsm}",
                f"!Sample_title = {title}",
                f"!Sample_description = {description}",
                f"!Sample_source_name_ch1 = {source}",
                f"!Sample_characteristics_ch1 = tissue: {source}",
                f"!Sample_relation = BioSample: https://example.invalid/{gsm}",
                f"!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term={srx}",
            ]
        )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def _write_tpm(path: Path, *, negative: bool = False) -> None:
    rows = [
        ("ENSG1|1|1|2|+|2|2|A|protein_coding", [1, 3, 0, 7]),
        ("ENSG2|1|3|4|+|2|2|A|protein_coding", [2, 5, 1, 1]),
        ("ENSG3|1|5|6|+|2|2|B|protein_coding", [4, 4, 8, 2]),
    ]
    if negative:
        rows[0][1][0] = -1
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["GeneID", "WM156A_OV", "WM156C_OM", "WM1859A_OV", "WM1859B_OM"])
        for gene_id, values in rows:
            writer.writerow([gene_id, *values])


def _source_manifest(path: Path, soft: Path, tpm: Path) -> None:
    data = {
        "schema_version": "hgsoc_external_sources.v1",
        "study_accession": "GSE277107",
        "files": [
            {
                "role": "geo_family_soft",
                "filename": soft.name,
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/test/soft.gz",
                "bytes": soft.stat().st_size,
                "sha256": sha256(soft),
            },
            {
                "role": "processed_tpm",
                "filename": tpm.name,
                "url": "https://ftp.ncbi.nlm.nih.gov/geo/test/tpm.gz",
                "bytes": tpm.stat().st_size,
                "sha256": sha256(tpm),
            },
        ],
        "expected": {
            "rna_samples": 4,
            "matched_pairs": 2,
            "source_gene_rows": 3,
            "unique_gene_symbols": 2,
        },
    }
    path.write_text(json.dumps(data), encoding="utf-8")


def test_pair_gate_preserves_submitted_anomaly(tmp_path: Path) -> None:
    soft = tmp_path / "family.soft.gz"
    _write_soft(soft)
    rows, anomalies = build_paired_metadata(soft, expected_samples=4, expected_pairs=2)
    assert {row["pair_id"] for row in rows} == {"WM156", "WM1859"}
    assert [row["normalized_site"] for row in rows if row["pair_id"] == "WM156"] == [
        "ovary",
        "omentum",
    ]
    assert anomalies == [
        {
            "geo_sample_accession": "GSM4",
            "field": "Sample_title",
            "value": "1859B_OM_HGSCS_omentum",
        }
    ]


def test_pair_gate_rejects_incomplete_pairs(tmp_path: Path) -> None:
    soft = tmp_path / "family.soft.gz"
    _write_soft(soft, omit_omentum=True)
    with pytest.raises(GSE277107Error, match="expected 4 GEO samples"):
        build_paired_metadata(soft, expected_samples=4, expected_pairs=2)


def test_prepare_aggregates_symbols_and_builds_paired_delta(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    soft = raw / "family.soft.gz"
    tpm = raw / "tpm.txt.gz"
    source_manifest = tmp_path / "sources.json"
    _write_soft(soft)
    _write_tpm(tpm)
    _source_manifest(source_manifest, soft, tpm)
    output = tmp_path / "processed"
    receipt = prepare_dataset(
        source_manifest_path=source_manifest, source_dir=raw, output_dir=output
    )
    assert receipt["status"] == "completed"
    assert receipt["validated_dimensions"] == {
        "rna_samples": 4,
        "matched_pairs": 2,
        "source_gene_rows": 3,
        "unique_gene_symbols": 2,
        "ovary_samples": 2,
        "omentum_samples": 2,
    }
    with gzip.open(output / "gene_symbol_tpm.tsv.gz", "rt", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0] == [
        "gene_name",
        "WM156A_OV",
        "WM156C_OM",
        "WM1859A_OV",
        "WM1859B_OM",
    ]
    assert [float(value) for value in rows[1][1:]] == [3.0, 8.0, 1.0, 8.0]
    with gzip.open(
        output / "paired_log2_tpm_delta_omentum_minus_ovary.tsv.gz",
        "rt",
        encoding="utf-8",
    ) as handle:
        delta = list(csv.reader(handle, delimiter="\t"))
    assert delta[0] == ["gene_name", "WM156", "WM1859"]
    assert float(delta[1][1]) == pytest.approx(math.log2(9) - math.log2(4), rel=1e-10)
    assert float(delta[1][2]) == pytest.approx(math.log2(9) - math.log2(2), rel=1e-10)
    loaded_manifest = REGULATORY_LOADER._load_manifest(
        output / "paired_sample_manifest.tsv", "GSE277107", False
    )
    loaded_samples, loaded_values = REGULATORY_LOADER._load_expression(
        output / "gene_symbol_tpm.tsv.gz"
    )
    assert {row["run_accession"] for row in loaded_manifest} == set(loaded_samples)
    assert loaded_values["A"][0] == pytest.approx(math.log1p(3.0))
    gate = audit_prepared_dataset(output)
    assert gate["scientific_success"] is True
    assert gate["regulatory_runner_contract"]["run_accession_matches_expression_columns"]


def test_receipt_gate_rejects_output_drift(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    soft = raw / "family.soft.gz"
    tpm = raw / "tpm.txt.gz"
    source_manifest = tmp_path / "sources.json"
    _write_soft(soft)
    _write_tpm(tpm)
    _source_manifest(source_manifest, soft, tpm)
    output = tmp_path / "processed"
    prepare_dataset(source_manifest_path=source_manifest, source_dir=raw, output_dir=output)
    with (output / "paired_sample_manifest.tsv").open("a", encoding="utf-8") as handle:
        handle.write("corruption\n")
    with pytest.raises(GSE277107Error, match="provenance mismatch"):
        audit_prepared_dataset(output)


def test_prepare_rejects_negative_tpm(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    soft = raw / "family.soft.gz"
    tpm = raw / "tpm.txt.gz"
    source_manifest = tmp_path / "sources.json"
    _write_soft(soft)
    _write_tpm(tpm, negative=True)
    _source_manifest(source_manifest, soft, tpm)
    with pytest.raises(GSE277107Error, match="finite and non-negative"):
        prepare_dataset(
            source_manifest_path=source_manifest,
            source_dir=raw,
            output_dir=tmp_path / "processed",
        )


def test_roihu_sbatch_is_syntax_valid_and_solver_free() -> None:
    script = ROOT / "hpc/roihu/gse277107_fetch_prepare.sbatch"
    subprocess.run(["bash", "-n", str(script)], check=True)
    text = script.read_text(encoding="utf-8")
    assert "audit_gse277107_receipt.py" in text
    assert "EXTERNAL_ROOT" in text
    assert "Gurobi" not in text
    assert "gurobi" not in text
