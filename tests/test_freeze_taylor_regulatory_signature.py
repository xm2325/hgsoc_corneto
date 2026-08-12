import argparse
import csv
import json
from pathlib import Path

from scripts.freeze_taylor_regulatory_signature import EXPECTED_STUDIES, freeze


def _edge(source: str, target: str) -> dict[str, object]:
    return {"source": source, "target": target, "sign": 1}


def _receipt(
    conditions: list[dict[str, object]], mode: str, study: str, bundle: str
) -> dict[str, object]:
    return {
        "status": "completed",
        "response_blind": True,
        "analysis_mode": mode,
        "study_accession": study,
        "bundle": {"sha256": bundle},
        "method": {
            "name": "CarnivalFlow",
            "single_joint_problem": True,
            "lambda_scaling": "mean_fit",
            "lambda_nominal": 0.001,
        },
        "scope_counts": {
            "included_conditions": len(conditions),
            "preprocessing_blocked": 0,
        },
        "conditions": conditions,
    }


def _write(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_freeze_uses_patient_prevalence_and_cohort_recurrence(tmp_path: Path) -> None:
    conditions = []
    index = 0
    cohort_conditions: dict[str, list[dict[str, object]]] = {}
    for study, count in EXPECTED_STUDIES.items():
        rows = []
        for _ in range(count):
            edges = []
            # Six patients across two cohorts: this edge passes 10% and 2/4.
            if index in {0, 1, 2, 9, 10, 11}:
                edges.append(_edge("A", "B"))
            # Five patients only: this edge fails the patient threshold.
            if index in {0, 1, 2, 9, 10}:
                edges.append(_edge("C", "D"))
            # Six patients, one cohort only: this edge fails cohort recurrence.
            if index in {3, 4, 5, 6, 7, 8}:
                edges.append(_edge("E", "F"))
            row = {
                "run_accession": f"R{index}",
                "patient_id": f"P{index % 52}",
                "study_accession": study,
                "status": "optimal",
                "selected_edges": edges,
            }
            rows.append(row)
            conditions.append(row)
            index += 1
        cohort_conditions[study] = rows

    pooled_path = tmp_path / "pooled.json"
    _write(pooled_path, _receipt(conditions, "pooled", "pooled", "bundle"))
    cohort_specs = []
    for study, rows in cohort_conditions.items():
        path = tmp_path / f"{study}.json"
        _write(path, _receipt(rows, "cohort", study, "bundle"))
        cohort_specs.append(f"{study}={path}")
    balanced_path = tmp_path / "balanced.json"
    _write(
        balanced_path,
        _receipt(conditions[:52], "pooled", "patient_balanced", "balanced"),
    )
    richer_path = tmp_path / "richer.json"
    _write(richer_path, _receipt(conditions, "pooled", "richer", "richer"))

    output = tmp_path / "signature.tsv"
    receipt_path = tmp_path / "receipt.json"
    result = freeze(
        argparse.Namespace(
            pooled=pooled_path,
            cohort=cohort_specs,
            patient_balanced=balanced_path,
            richer_pooled=richer_path,
            min_patient_fraction=0.10,
            min_cohort_fraction=0.50,
            output=output,
            receipt=receipt_path,
        )
    )
    assert result["selected_edges"] == 1
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["feature_id"] == "A|B|+"
    assert rows[0]["taylor_patient_count"] == "6"
    assert rows[0]["cohort_recurrence"] == "2"
