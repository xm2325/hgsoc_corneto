import argparse
import csv
import hashlib
import json
from pathlib import Path

from scripts.build_external_regulatory_evidence import build


def _tsv(path: Path, rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_build_patient_level_complete_frozen_grid(tmp_path: Path) -> None:
    signature = tmp_path / "signature.tsv"
    _tsv(
        signature,
        [
            ["feature_type", "feature_id", "expected_direction", "source", "target", "sign"],
            ["edge", "A|B|+", 1, "A", "B", 1],
        ],
    )
    manifest = tmp_path / "manifest.tsv"
    _tsv(
        manifest,
        [
            ["study_accession", "run_accession", "patient", "site"],
            ["TEST", "R1", "P1", "ovary"],
            ["TEST", "R2", "P2", "ovary"],
            ["TEST", "R3", "P1", "omentum"],
            ["TEST", "R4", "P2", "omentum"],
        ],
    )
    bundle = tmp_path / "bundle.json"
    bundle.write_text(
        json.dumps(
            {
                "sources": {"required_signature": {"sha256": _sha(signature)}},
                "graph": [{"source": "A", "target": "B", "sign": 1}],
            }
        ),
        encoding="utf-8",
    )
    edge = {"source": "A", "target": "B", "sign": 1}
    solution = tmp_path / "solution.json"
    solution.write_text(
        json.dumps(
            {
                "status": "completed",
                "response_blind": True,
                "method": {"lambda_nominal": 0.001, "lambda_scaling": "mean_fit"},
                "solver": {"selected": "GUROBI", "has_incumbent": True},
                "bundle": {"sha256": _sha(bundle)},
                "conditions": [
                    {
                        "run_accession": run,
                        "status": "optimal",
                        "selected_edges": [edge] if run in {"R1", "R3"} else [],
                    }
                    for run in ("R1", "R2", "R3", "R4")
                ],
            }
        ),
        encoding="utf-8",
    )
    normalization = tmp_path / "normalization.json"
    normalization.write_text('{"status":"completed"}\n', encoding="utf-8")
    output = tmp_path / "evidence"
    result = build(
        argparse.Namespace(
            solution=solution,
            bundle=bundle,
            signature=signature,
            manifest=manifest,
            normalization_receipt=normalization,
            study="TEST",
            patient_id_field="patient",
            analysis_unit="patient_tissue",
            group=["ovary=site:ovary", "omentum=site:omentum"],
            within_patient_threshold=0.5,
            output_dir=output,
        )
    )
    assert result["groups"]["ovary"]["patients"] == 2
    with (output / "ovary.evidence.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert {row["evaluable"] for row in rows} == {"1"}
    assert json.loads((output / "omentum.contract.json").read_text())[
        "inference"
    ]["all_frozen_edges_present_in_candidate_graph"] is True
