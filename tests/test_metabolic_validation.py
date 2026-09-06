from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from hgsoc_corneto.metabolic.validation import (
    audit_primal,
    optimization_issues,
    require_independent_receipts,
    sha256,
    validate_receipt,
)


def good_telemetry():
    return {
        "cvxpy_status": "optimal",
        "gurobi_status_code": 2,
        "gurobi_status_name": "OPTIMAL",
        "has_incumbent": True,
        "solution_count": 1,
        "objective_value": -10.0,
        "best_bound": -10.0,
        "relative_gap": 0.0,
        "requested_mip_gap": 1e-4,
    }


def instrument_receipt(tmp_path, receipt, summaries):
    """Make a small complete optimization certificate for gate tests."""
    artifacts = {}
    for label in ("solution", "mip_start", "gurobi_log"):
        path = tmp_path / label
        path.write_text("test solver artifact\n")
        artifacts[label] = str(path)
    receipt["scientific_success"] = True
    telemetry = good_telemetry()
    telemetry.update(
        objective_value=summaries[0]["problem_objective_value"],
        best_bound=summaries[0]["problem_objective_value"],
    )
    receipt["instrumentation"] = {
        "telemetry": telemetry,
        "summary_error": None,
        "summaries": summaries,
        "artifacts": artifacts,
        "artifact_sha256": {k: sha256(Path(v)) for k, v in artifacts.items()},
        "primal_validation": {
            "status": "passed",
            "conditions": [
                {
                    "condition": s["condition"],
                    "passed": True,
                    "max_mass_balance_residual": 0.0,
                    "max_bound_violation": 0.0,
                    "max_integrality_error": 0.0,
                    "flux_on_unselected_indicator_count": 0,
                }
                for s in summaries
            ],
        },
    }
    return receipt


def good_receipt(tmp_path):
    summary = {
        "condition": "A",
        "status": "optimal",
        "problem_objective_value": -10.0,
        "active_by_flux": ["R"],
        "active_by_indicator": ["R"],
        "nonzero_fluxes": [["R", 1.0]],
    }
    return instrument_receipt(
        tmp_path,
        {
            "status": "completed",
            "schema_version": "metabolic_independent_checkpoint.v1",
            "context_sha256": "context",
            "condition": "A",
            "solution": summary,
        },
        [summary],
    )


def test_optimal_inaccurate_and_time_limit_do_not_pass():
    telemetry = good_telemetry()
    assert optimization_issues(telemetry) == []
    telemetry.update(
        cvxpy_status="optimal_inaccurate",
        gurobi_status_code=9,
        gurobi_status_name="TIME_LIMIT",
        relative_gap=0.035,
        best_bound=-10.35,
    )
    assert "gap_exceeds_requested_tolerance" in optimization_issues(telemetry)
    assert "cvxpy_status_not_optimal" in optimization_issues(telemetry)


@pytest.mark.parametrize(
    "key,value",
    [
        ("relative_gap", float("nan")),
        ("best_bound", float("inf")),
        ("objective_value", None),
        ("solution_count", 0),
    ],
)
def test_incomplete_telemetry_fails_closed(key, value):
    telemetry = good_telemetry()
    telemetry[key] = value
    assert optimization_issues(telemetry)


def test_recomputes_gap_from_objective_and_bound():
    telemetry = good_telemetry()
    telemetry["best_bound"] = -11.0  # reported relative_gap=0 cannot hide this
    assert "objective_bound_gap_exceeds_requested_tolerance" in optimization_issues(telemetry)


@pytest.mark.parametrize("fault", ["gap", "context", "primal", "artifact", "nan", "indicator"])
def test_rejects_incorrect_receipts(tmp_path, fault):
    receipt = good_receipt(tmp_path)
    if fault == "gap":
        receipt["instrumentation"]["telemetry"]["relative_gap"] = 0.03
    elif fault == "context":
        receipt["context_sha256"] = "other"
    elif fault == "primal":
        del receipt["instrumentation"]["primal_validation"]
    elif fault == "artifact":
        (tmp_path / "solution").write_text("modified after receipt was written")
    elif fault == "nan":
        receipt["solution"]["nonzero_fluxes"][0][1] = float("nan")
    elif fault == "indicator":
        receipt["solution"]["active_by_indicator"] = []
    path = tmp_path / "receipt.json"
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="receipt rejected"):
        validate_receipt(path, context_sha256="context", conditions=["A"], kind="independent")


def test_valid_receipt_and_missing_cohort_member(tmp_path):
    context = tmp_path / "context.json"
    context.write_text(json.dumps({"conditions": ["A", "B"]}))
    directory = tmp_path / "independent"
    directory.mkdir()
    receipt = good_receipt(tmp_path)
    receipt["context_sha256"] = sha256(context)
    path = directory / "000_A.json"
    path.write_text(json.dumps(receipt))
    validate_receipt(path, context_sha256=sha256(context), conditions=["A"], kind="independent")
    with pytest.raises(FileNotFoundError):
        require_independent_receipts(context)


class Metabolite:
    id = "A"


def primal_model(flow=(1.0, 1.0), indicators=(1.0, 1.0)):
    m = Metabolite()
    reactions = {
        "in": SimpleNamespace(lower_bound=0.0, upper_bound=1000.0, metabolites={m: 1.0}),
        "out": SimpleNamespace(lower_bound=0.0, upper_bound=1000.0, metabolites={m: -1.0}),
    }
    model = SimpleNamespace(reactions=SimpleNamespace(get_by_id=reactions.__getitem__))
    problem = SimpleNamespace(
        expr=SimpleNamespace(
            flow=SimpleNamespace(value=np.asarray(flow)),
            edge_has_flux=SimpleNamespace(value=np.asarray(indicators)),
        )
    )
    return model, problem


def test_trickle_flow_fails_even_when_raw_mass_balance_passes():
    model, problem = primal_model((0.009, 0.009), (9e-6, 9e-6))
    result = audit_primal(model, problem, ["in", "out"], ["A"], {"A": {}})
    assert result["status"] == "failed"
    assert result["conditions"][0]["max_mass_balance_residual"] == 0
    assert result["conditions"][0]["flux_on_unselected_indicator_count"] == 2


def test_primal_checks_full_flux_and_context_bounds():
    model, problem = primal_model()
    assert audit_primal(model, problem, ["in", "out"], ["A"], {"A": {}})["status"] == "passed"
    result = audit_primal(model, problem, ["in", "out"], ["A"], {"A": {"out": (0, 0.5)}})
    assert result["status"] == "failed"
    assert result["conditions"][0]["max_bound_violation"] == 0.5


def test_joint_checks_all_independent_receipts_before_importing_solver(tmp_path):
    script = Path(__file__).resolve().parents[1] / "scripts/run_corneto_metabolic_instrumented.py"
    spec = importlib.util.spec_from_file_location("instrumented_runner", script)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    context = tmp_path / "context.json"
    context.write_text(json.dumps({"status": "prepared", "conditions": ["A"]}))
    args = SimpleNamespace(context=context, output=tmp_path / "joint.json")
    with pytest.raises(FileNotFoundError):
        runner.joint(args)


def test_scientific_review_hold_stops_before_solver_import(tmp_path):
    script = Path(__file__).resolve().parents[1] / "scripts/run_corneto_metabolic_instrumented.py"
    spec = importlib.util.spec_from_file_location("held_runner", script)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    context = tmp_path / "context.json"
    context.write_text(json.dumps({"status": "prepared", "conditions": ["A"]}))
    hold = tmp_path / "scientific_review_hold.json"
    hold.write_text(json.dumps({"status": "hold", "reason": "input review"}))
    with pytest.raises(RuntimeError, match="solver was not started"):
        runner._context(context)
    hold.write_text(json.dumps({"status": "released"}))
    assert runner._context(context)[0]["status"] == "prepared"
