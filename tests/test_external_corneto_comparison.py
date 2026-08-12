import csv
import hashlib
import json
from pathlib import Path

import pytest

from scripts.compare_external_corneto import ExternalValidationError, compare


def _tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path: Path) -> tuple[Path, list[str], dict[str, Path]]:
    signature = tmp_path / "signature.tsv"
    _tsv(
        signature,
        [
            {"feature_type": "edge", "feature_id": "A|B|+", "expected_direction": 1},
            {"feature_type": "reaction", "feature_id": "MAR001", "expected_direction": -1},
        ],
    )
    groups = {
        "tumour": {
            "p1": [(True, 1), (True, -1)],
            "p2": [(True, 1), (False, 0)],
            "p3": [(True, -1), (True, -1)],
        },
        "normal": {
            "p1": [(True, 1), (False, 0)],
            "p4": [(False, 0), (False, 0)],
            "p5": [(False, 0), (True, 1)],
        },
    }
    paths: dict[str, Path] = {}
    specs = []
    features = [("edge", "A|B|+"), ("reaction", "MAR001")]
    for label, patients in groups.items():
        evidence = tmp_path / f"{label}.tsv"
        rows = []
        for patient, values in patients.items():
            for (feature_type, feature_id), (selected, direction) in zip(
                features, values, strict=True
            ):
                rows.append(
                    {
                        "patient_id": patient,
                        "feature_type": feature_type,
                        "feature_id": feature_id,
                        "selected": int(selected),
                        "direction": direction,
                    }
                )
        _tsv(evidence, rows)
        contract = tmp_path / f"{label}.json"
        contract.write_text(
            json.dumps(
                {
                    "schema_version": "external_corneto_group.v1",
                    "status": "completed",
                    "group_label": label,
                    "source_accession": "GSE-test",
                    "evidence_sha256": _sha(evidence),
                    "signature_sha256": _sha(signature),
                    "patient_count": len(patients),
                    "analysis_unit": "patient_pseudobulk",
                    "normalization": {
                        "performed_within_dataset": True,
                        "pooled_raw_expression": False,
                        "input_scale": "network_selection",
                    },
                    "independence": {
                        "cells_as_replicates": False,
                        "patient_id_column": "patient_id",
                    },
                    "inference": {
                        "signature_frozen_before_external_scoring": True,
                        "feature_selection_using_external_labels": False,
                    },
                }
            ),
            encoding="utf-8",
        )
        paths[label] = contract
        specs.append(f"{label}={evidence},{contract}")
    return signature, specs, paths


def test_patient_level_external_comparison(tmp_path: Path) -> None:
    signature, specs, _ = _fixture(tmp_path)
    result = compare(
        signature_path=signature,
        group_specs=specs,
        bootstrap_iterations=200,
        seed=7,
    )
    assert result["status"] == "completed"
    assert result["contract"]["bootstrap_unit"] == "patient_id"
    assert result["groups"]["tumour"]["patient_count"] == 3
    features = result["groups"]["tumour"]["features"]
    assert features[0]["prevalence"] == 1.0
    assert features[0]["direction_concordance"] == pytest.approx(2 / 3)
    pair = result["pairwise"][0]
    assert pair["consensus_jaccard"] == 0.0
    assert pair["matched_patient_count"] == 1
    assert pair["mean_matched_patient_jaccard"] == 0.5


def test_raw_expression_pooling_fails_closed(tmp_path: Path) -> None:
    signature, specs, contracts = _fixture(tmp_path)
    contract = json.loads(contracts["normal"].read_text(encoding="utf-8"))
    contract["normalization"]["pooled_raw_expression"] = True
    contracts["normal"].write_text(json.dumps(contract), encoding="utf-8")
    with pytest.raises(ExternalValidationError, match="pooled raw expression"):
        compare(signature_path=signature, group_specs=specs, bootstrap_iterations=100)


def test_cells_as_n_fails_closed(tmp_path: Path) -> None:
    signature, specs, contracts = _fixture(tmp_path)
    contract = json.loads(contracts["tumour"].read_text(encoding="utf-8"))
    contract["independence"]["cells_as_replicates"] = True
    contracts["tumour"].write_text(json.dumps(contract), encoding="utf-8")
    with pytest.raises(ExternalValidationError, match="cells-as-n"):
        compare(signature_path=signature, group_specs=specs, bootstrap_iterations=100)


def test_incomplete_patient_signature_grid_fails_closed(tmp_path: Path) -> None:
    signature, specs, contracts = _fixture(tmp_path)
    evidence = Path(specs[0].split("=", 1)[1].split(",", 1)[0])
    rows = list(csv.DictReader(evidence.open(encoding="utf-8"), delimiter="\t"))
    _tsv(evidence, rows[:-1])
    contract = json.loads(contracts["tumour"].read_text(encoding="utf-8"))
    contract["evidence_sha256"] = _sha(evidence)
    contracts["tumour"].write_text(json.dumps(contract), encoding="utf-8")
    with pytest.raises(ExternalValidationError, match="complete frozen grid"):
        compare(signature_path=signature, group_specs=specs, bootstrap_iterations=100)
