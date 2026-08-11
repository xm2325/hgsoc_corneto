import csv
import json
from pathlib import Path

import pytest

from hgsoc_corneto.phenotype import (
    PhenotypeIntakeError,
    blocked_receipt,
    validate_phenotype_intake,
)

FIELDS = [
    "canonical_ocm_id",
    "patient_id",
    "drug",
    "endpoint",
    "value",
    "unit",
    "endpoint_definition",
    "is_exact",
    "source_file",
    "source_record_id",
]


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _inputs(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    _write_tsv(
        manifest,
        ["run_accession", "canonical_ocm_id", "patient_id", "primary_cohort_eligible"],
        [
            {
                "run_accession": "R1",
                "canonical_ocm_id": "OCM1",
                "patient_id": "OCM1",
                "primary_cohort_eligible": "true",
            },
            {
                "run_accession": "R2",
                "canonical_ocm_id": "OCM2-1",
                "patient_id": "OCM2",
                "primary_cohort_eligible": "true",
            },
        ],
    )
    schema = tmp_path / "schema.json"
    schema.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "cohort": {"primary_ocms": 2},
                "endpoints": {
                    "paclitaxel_auc": {
                        "role": "primary",
                        "drug": "paclitaxel",
                        "constraint": "nonnegative",
                        "accepted_units": ["source_reported_auc"],
                        "transform_after_gate": None,
                    },
                    "paclitaxel_gi50": {
                        "role": "secondary",
                        "drug": "paclitaxel",
                        "constraint": "positive",
                        "accepted_units": ["nM"],
                        "transform_after_gate": "log10",
                    },
                    "cumulative_paclitaxel_exposure": {
                        "role": "exposure",
                        "drug": "paclitaxel",
                        "constraint": "nonnegative",
                        "accepted_units": ["mg"],
                        "transform_after_gate": None,
                    },
                },
                "analysis_readiness": {
                    "primary_endpoint": "paclitaxel_auc",
                    "secondary_endpoint": "paclitaxel_gi50",
                    "exposure_endpoint": "cumulative_paclitaxel_exposure",
                },
            }
        ),
        encoding="utf-8",
    )
    return manifest, schema


def _row(ocm: str, patient: str, endpoint: str, value: str, unit: str) -> dict[str, str]:
    return {
        "canonical_ocm_id": ocm,
        "patient_id": patient,
        "drug": "paclitaxel",
        "endpoint": endpoint,
        "value": value,
        "unit": unit,
        "endpoint_definition": f"exact definition for {endpoint}",
        "is_exact": "true",
        "source_file": "Taylor_underlying.xlsx",
        "source_record_id": f"{ocm}:{endpoint}",
    }


def test_missing_file_is_fail_closed(tmp_path):
    receipt = blocked_receipt(
        phenotype_path=tmp_path / "missing.tsv",
        manifest_path=tmp_path / "manifest.tsv",
        schema_path=tmp_path / "schema.json",
    )
    assert receipt["status"] == "blocked_missing_phenotype_file"
    assert receipt["association_allowed"] is False
    assert receipt["association_run"] is False


def test_complete_exact_table_is_ready_without_values_in_receipt(tmp_path):
    manifest, schema = _inputs(tmp_path)
    rows = []
    for ocm, patient in (("OCM1", "OCM1"), ("OCM2-1", "OCM2")):
        rows.extend(
            [
                _row(ocm, patient, "paclitaxel_auc", "0.4", "source_reported_auc"),
                _row(ocm, patient, "paclitaxel_gi50", "5.0", "nM"),
                _row(ocm, patient, "cumulative_paclitaxel_exposure", "0", "mg"),
            ]
        )
    phenotype = tmp_path / "phenotype.tsv"
    _write_tsv(phenotype, FIELDS, rows)
    receipt = validate_phenotype_intake(
        phenotype_path=phenotype, manifest_path=manifest, schema_path=schema
    )
    assert receipt["status"] == "ready"
    assert receipt["association_allowed"] is True
    assert receipt["endpoints"]["paclitaxel_gi50"]["transform_after_gate"] == "log10"
    assert "0.4" not in json.dumps(receipt)
    assert receipt["safety"]["association_run"] is False


def test_incomplete_primary_endpoint_is_blocked(tmp_path):
    manifest, schema = _inputs(tmp_path)
    phenotype = tmp_path / "phenotype.tsv"
    _write_tsv(
        phenotype,
        FIELDS,
        [_row("OCM1", "OCM1", "paclitaxel_auc", "0.4", "source_reported_auc")],
    )
    receipt = validate_phenotype_intake(
        phenotype_path=phenotype, manifest_path=manifest, schema_path=schema
    )
    assert receipt["status"] == "blocked_incomplete_exact_phenotype"
    assert receipt["association_allowed"] is False
    assert receipt["endpoints"]["paclitaxel_auc"]["missing_ocm_ids"] == ["OCM2-1"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"patient_id": "WRONG"}, "patient mapping conflict"),
        ({"unit": "uM"}, "unit 'uM' is not accepted"),
        ({"is_exact": "false"}, "proxy/imputed rows are forbidden"),
    ],
)
def test_unsafe_rows_are_rejected(tmp_path, mutation, message):
    manifest, schema = _inputs(tmp_path)
    row = _row("OCM1", "OCM1", "paclitaxel_gi50", "5", "nM")
    row.update(mutation)
    phenotype = tmp_path / "phenotype.tsv"
    _write_tsv(phenotype, FIELDS, [row])
    with pytest.raises(PhenotypeIntakeError, match=message):
        validate_phenotype_intake(
            phenotype_path=phenotype, manifest_path=manifest, schema_path=schema
        )


def test_duplicate_ocm_endpoint_is_rejected(tmp_path):
    manifest, schema = _inputs(tmp_path)
    row = _row("OCM1", "OCM1", "paclitaxel_auc", "0.4", "source_reported_auc")
    duplicate = {**row, "source_record_id": "second"}
    phenotype = tmp_path / "phenotype.tsv"
    _write_tsv(phenotype, FIELDS, [row, duplicate])
    with pytest.raises(PhenotypeIntakeError, match="duplicate OCM/endpoint"):
        validate_phenotype_intake(
            phenotype_path=phenotype, manifest_path=manifest, schema_path=schema
        )
