import csv
import gzip
import hashlib

import pytest

from hgsoc_corneto.external_validation import (
    audit_gse_count_matrix,
    extract_depmap_gene_effects,
    extract_gse_candidate_log_cpm,
    normalize_cell_line_name,
    parse_depmap_gene_label,
    resolve_depmap_hgsoc_models,
    summarize_dependency_rows,
    summarize_expression_rows,
)


def test_audit_gse_matrix_is_checksum_and_schema_gated(tmp_path):
    path = tmp_path / "counts.tsv.gz"
    content = "\tPDO1\tFT1\nENSG1\t1\t0\nENSG2\t3\t4\n"
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as handle:
            handle.write(content.encode())
    checksum = hashlib.sha256(path.read_bytes()).hexdigest()
    result = audit_gse_count_matrix(
        path,
        expected_sha256=checksum,
        expected_samples=["PDO1", "FT1"],
        expected_gene_rows=2,
    )
    assert result["scientific_success"] is True
    assert result["gene_row_count"] == 2


def test_audit_gse_matrix_rejects_unexpected_sample(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text("\tPDO1\tEXTRA\nENSG1\t1\t0\n", encoding="utf-8")
    checksum = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="sample set mismatch"):
        audit_gse_count_matrix(
            path,
            expected_sha256=checksum,
            expected_samples=["PDO1", "FT1"],
            expected_gene_rows=1,
        )


def test_resolve_depmap_models_never_promotes_all_ovarian_lines():
    models = [
        {
            "ModelID": "ACH-CORE",
            "StrippedCellLineName": "KURAMOCHI",
            "OncotreeLineage": "Ovary/Fallopian Tube",
        },
        {
            "ModelID": "ACH-OTHER",
            "StrippedCellLineName": "SKOV3",
            "OncotreeLineage": "Ovary/Fallopian Tube",
        },
        {"ModelID": "ACH-LUNG", "StrippedCellLineName": "A549", "OncotreeLineage": "Lung"},
    ]
    curated = [
        {
            "model_id": "ACH-CORE",
            "cell_line": "KURAMOCHI",
            "evidence_tier": "high_confidence_hgsoc_like",
            "evidence_url": "https://example.test/evidence",
        }
    ]
    positives, comparators = resolve_depmap_hgsoc_models(models, curated)
    assert [row["ModelID"] for row in positives] == ["ACH-CORE"]
    assert [row["ModelID"] for row in comparators] == ["ACH-OTHER"]
    assert comparators[0]["validation_group"] == "other_ovarian_not_hgsoc_positive"


def test_resolve_depmap_models_fails_on_release_identity_drift():
    models = [
        {
            "ModelID": "ACH-CORE",
            "StrippedCellLineName": "A_DIFFERENT_LINE",
            "OncotreeLineage": "Ovary/Fallopian Tube",
        }
    ]
    curated = [
        {
            "model_id": "ACH-CORE",
            "cell_line": "KURAMOCHI",
            "evidence_tier": "high_confidence_hgsoc_like",
            "evidence_url": "https://example.test/evidence",
        }
    ]
    with pytest.raises(ValueError, match="identity mismatch"):
        resolve_depmap_hgsoc_models(models, curated)


def test_extract_and_summarize_depmap_effects(tmp_path):
    matrix = tmp_path / "CRISPRGeneEffect.csv"
    with matrix.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ModelID", "TPI1 (7167)", "TP53 (7157)"])
        writer.writerow(["ACH-CORE", "-1.2", "-0.1"])
        writer.writerow(["ACH-OTHER", "-0.2", "-0.3"])
    rows, missing_models, missing_genes = extract_depmap_gene_effects(
        matrix,
        model_groups={
            "ACH-CORE": "hgsoc_like_positive",
            "ACH-OTHER": "other_ovarian_not_hgsoc_positive",
        },
        candidate_symbols=["TPI1"],
    )
    assert not missing_models
    assert not missing_genes
    assert [row["gene_effect"] for row in rows] == [-1.2, -0.2]
    summary = summarize_dependency_rows(rows)
    assert {row["validation_group"] for row in summary} == {
        "hgsoc_like_positive",
        "other_ovarian_not_hgsoc_positive",
    }


def test_depmap_label_and_model_name_normalization():
    assert parse_depmap_gene_label("TPI1 (7167)") == "TPI1"
    assert normalize_cell_line_name("OVCAR-4_OVARY") == "OVCAR4"


def test_extract_gse_candidate_expression_uses_explicit_mapping(tmp_path):
    matrix = tmp_path / "counts.tsv"
    matrix.write_text(
        "\tPDO1\tFT1\nENSG_TPI1\t9\t1\nENSG_OTHER\t91\t99\n",
        encoding="utf-8",
    )
    rows, missing = extract_gse_candidate_log_cpm(
        matrix,
        sample_groups={"PDO1": "hgsoc_organoid", "FT1": "fallopian_tube_organoid"},
        gene_id_to_symbol={"ENSG_TPI1": "TPI1"},
        candidate_symbols=["TPI1"],
    )
    assert not missing
    assert [row["raw_count"] for row in rows] == [9, 1]
    summary = summarize_expression_rows(rows)
    assert {row["validation_group"] for row in summary} == {
        "hgsoc_organoid",
        "fallopian_tube_organoid",
    }


def test_extract_gse_candidate_expression_preserves_missing_gene(tmp_path):
    matrix = tmp_path / "counts.tsv"
    matrix.write_text("\tPDO1\nENSG_OTHER\t1\n", encoding="utf-8")
    rows, missing = extract_gse_candidate_log_cpm(
        matrix,
        sample_groups={"PDO1": "hgsoc_organoid"},
        gene_id_to_symbol={"ENSG_OTHER": "OTHER"},
        candidate_symbols=["TPI1"],
    )
    assert rows == []
    assert missing == ["TPI1"]
