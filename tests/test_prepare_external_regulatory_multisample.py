import argparse
import csv
from pathlib import Path

from scripts.prepare_external_regulatory_multisample import build_bundle


def _tsv(path: Path, rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def test_external_bundle_is_independently_standardized_and_patient_aware(
    tmp_path: Path,
) -> None:
    expression = tmp_path / "expression.tsv"
    _tsv(
        expression,
        [
            ["gene_name", "S1", "S2", "S3"],
            ["GENE1", 10, 1, 5],
            ["SRC", 1, 10, 5],
        ],
    )
    manifest = tmp_path / "manifest.tsv"
    _tsv(
        manifest,
        [
            ["study_accession", "run_accession", "patient_id", "comparison_role"],
            ["GSETEST", "S1", "P1", "tumour"],
            ["GSETEST", "S2", "P2", "reference"],
            ["GSETEST", "S3", "P3", "reference"],
        ],
    )
    collectri = tmp_path / "collectri.tsv"
    _tsv(collectri, [["source", "target", "sign"], ["TF1", "GENE1", 1]])
    pkn = tmp_path / "pkn.tsv"
    _tsv(pkn, [["source", "target", "sign"], ["SRC", "TF1", 1]])
    signature = tmp_path / "signature.tsv"
    _tsv(
        signature,
        [
            ["feature_type", "feature_id", "source", "target", "sign"],
            ["edge", "X|Y|+", "X", "Y", 1],
        ],
    )
    output = tmp_path / "bundle.json"
    result = build_bundle(
        argparse.Namespace(
            expression=expression,
            manifest=manifest,
            collectri=collectri,
            pkn=pkn,
            required_signature=signature,
            study="GSETEST",
            output=output,
            expected_count=3,
            include_role=["tumour", "reference"],
            role_field="comparison_role",
            patient_id_field="patient_id",
            min_targets=1,
            max_inputs=1,
            max_outputs=1,
            max_depth=2,
            max_edges=10,
        )
    )
    assert result["status"] == "completed"
    assert result["input_counts"]["selected_conditions"] == 3
    assert result["input_counts"]["unique_patients"] == 3
    assert result["input_counts"]["frozen_graph_union_edges"] == 2
    assert result["input_counts"]["required_signature_edges"] == 1
    assert {row["patient_id"] for row in result["conditions"]} == {"P1", "P2", "P3"}
