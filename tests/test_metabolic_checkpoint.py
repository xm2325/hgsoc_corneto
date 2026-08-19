from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_corneto_metabolic_checkpoint.py"


def _write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


def test_assemble_requires_complete_matching_checkpoints(tmp_path: Path) -> None:
    conditions = ["RUN1", "RUN2"]
    context = {
        "status": "prepared",
        "schema_version": "metabolic_checkpoint_context.v1",
        "study_accession": "TEST",
        "primary_only": True,
        "sample_count": 2,
        "samples": [{"run_accession": value} for value in conditions],
        "conditions": conditions,
        "expression": {"path": "matrix.tsv.gz", "transform": "log1p_tpm", "gene_count": 5},
        "model": {"path": "model.xml", "reactions": 4, "genes": 3, "biomass_id": "biomass_human"},
        "solver": {
            "requested": "gurobi",
            "used": "gurobi",
            "available": ["GUROBI"],
            "fallback_reason": None,
            "solve_fallback": None,
        },
        "candidate_selection": {"max_candidates": 25, "selected_count": 1},
        "objective": {"growth_fraction": 0.9, "independent_lambda": 0.1, "joint_lambda": 1.0},
        "repo_commit": "abc",
    }
    context_path = tmp_path / "context.json"
    _write(context_path, context)
    import hashlib

    context_sha = hashlib.sha256(context_path.read_bytes()).hexdigest()
    independent_dir = tmp_path / "independent"
    solutions = []
    for index, condition in enumerate(conditions):
        solution = {
            "condition": condition,
            "status": "optimal",
            "problem_objective_value": 1.0,
            "active_by_flux": [f"R{index}", "biomass_human"],
            "active_by_indicator": [f"R{index}", "biomass_human"],
            "nonzero_fluxes": [[f"R{index}", 1.0], ["biomass_human", 0.9]],
        }
        solutions.append(solution)
        _write(
            independent_dir / f"{index:03d}_{condition}.json",
            {
                "status": "completed",
                "schema_version": "metabolic_independent_checkpoint.v1",
                "condition": condition,
                "context_sha256": context_sha,
                "solution": solution,
                "slurm_job_id": str(100 + index),
            },
        )
    joint_path = tmp_path / "joint.json"
    _write(
        joint_path,
        {
            "status": "completed",
            "schema_version": "metabolic_joint_checkpoint.v1",
            "context_sha256": context_sha,
            "slurm_job_id": "200",
            "result": {
                "conditions": conditions,
                "active_tolerance": 1e-7,
                "joint": solutions,
            },
        },
    )
    output = tmp_path / "final.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "assemble",
            "--context",
            str(context_path),
            "--independent-dir",
            str(independent_dir),
            "--joint-receipt",
            str(joint_path),
            "--output",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["status"] == "completed"
    assert result["sample_count"] == 2
    assert result["corneto"]["conditions"] == conditions
    assert result["corneto"]["independent_active_union_size"] == 3
    assert result["corneto"]["joint_active_union_size"] == 3
    assert result["checkpoint_provenance"]["context_sha256"] == context_sha


def test_assemble_fails_closed_when_a_checkpoint_is_missing(tmp_path: Path) -> None:
    context_path = tmp_path / "context.json"
    _write(
        context_path,
        {
            "status": "prepared",
            "schema_version": "metabolic_checkpoint_context.v1",
            "study_accession": "TEST",
            "conditions": ["RUN1"],
        },
    )
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "assemble",
            "--context",
            str(context_path),
            "--independent-dir",
            str(tmp_path / "missing"),
            "--joint-receipt",
            str(tmp_path / "joint.json"),
            "--output",
            str(tmp_path / "must_not_exist.json"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0
    assert not (tmp_path / "must_not_exist.json").exists()
