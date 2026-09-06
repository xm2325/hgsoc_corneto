"""Fail-closed acceptance of instrumented optimization evidence.

Passing these checks certifies the recorded optimization/numerical contract,
not biological validity of expression constraints or culture-medium choices.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def optimization_issues(telemetry: Mapping[str, Any]) -> list[str]:
    issues = []
    if telemetry.get("cvxpy_status") != "optimal":
        issues.append("cvxpy_status_not_optimal")
    if telemetry.get("gurobi_status_code") != 2 or telemetry.get("gurobi_status_name") != "OPTIMAL":
        issues.append("gurobi_status_not_optimal")
    count = telemetry.get("solution_count")
    if telemetry.get("has_incumbent") is not True or not _finite(count) or count < 1:
        issues.append("missing_incumbent")
    for key in ("objective_value", "best_bound", "relative_gap", "requested_mip_gap"):
        if not _finite(telemetry.get(key)):
            issues.append(f"missing_or_nonfinite_{key}")
    gap, target = telemetry.get("relative_gap"), telemetry.get("requested_mip_gap")
    if _finite(gap) and _finite(target):
        if not 0 <= target < 1 or gap < 0 or gap > target + 1e-12:
            issues.append("gap_exceeds_requested_tolerance")
    objective, bound = telemetry.get("objective_value"), telemetry.get("best_bound")
    if _finite(objective) and _finite(bound) and _finite(target):
        actual = (
            abs(objective - bound) / abs(objective)
            if objective
            else (0.0 if bound == 0 else math.inf)
        )
        if actual > target + 1e-12:
            issues.append("objective_bound_gap_exceeds_requested_tolerance")
    return issues


def audit_primal(
    cobra_model: Any,
    problem: Any,
    reaction_ids: Sequence[str],
    conditions: Sequence[str],
    reaction_bounds: Mapping[str, Mapping[str, tuple[float | None, float | None]]],
    *,
    feasibility_tolerance: float = 1e-6,
    activity_tolerance: float = 1e-7,
) -> dict[str, Any]:
    """Check full-precision primal fluxes before summary truncation."""
    import numpy as np

    flows = np.asarray(problem.expr.flow.value, dtype=float)
    indicators = np.asarray(problem.expr.edge_has_flux.value, dtype=float)
    shape = (len(reaction_ids), len(conditions))
    if len(conditions) == 1:
        flows, indicators = flows.reshape(-1, 1), indicators.reshape(-1, 1)
    if flows.shape != shape or indicators.shape != shape:
        return {"status": "failed", "issues": ["invalid_primal_shape"]}
    if not np.isfinite(flows).all() or not np.isfinite(indicators).all():
        return {"status": "failed", "issues": ["nonfinite_primal"]}
    rows = []
    for column, condition in enumerate(conditions):
        residuals: dict[str, float] = {}
        bound_max = 0.0
        for index, rid in enumerate(reaction_ids):
            reaction = cobra_model.reactions.get_by_id(rid)
            lower, upper = float(reaction.lower_bound), float(reaction.upper_bound)
            override = reaction_bounds[condition].get(rid, (None, None))
            if override[0] is not None:
                lower = max(lower, float(override[0]))
            if override[1] is not None:
                upper = min(upper, float(override[1]))
            value = float(flows[index, column])
            bound_max = max(bound_max, lower - value, value - upper)
            for metabolite, coefficient in reaction.metabolites.items():
                residuals[metabolite.id] = residuals.get(metabolite.id, 0.0) + coefficient * value
        y = indicators[:, column]
        off_flux = np.abs(flows[:, column])[(y < 0.5)]
        mass_max = max(map(abs, residuals.values()), default=0.0)
        integer_error = float(np.max(np.abs(y - np.round(y)), initial=0))
        off_count = int(np.sum(off_flux > activity_tolerance))
        passed = (
            mass_max <= feasibility_tolerance
            and bound_max <= feasibility_tolerance
            and integer_error <= 1e-5
            and off_count == 0
            and float(np.min(y, initial=0)) >= -1e-5
            and float(np.max(y, initial=1)) <= 1 + 1e-5
        )
        rows.append(
            {
                "condition": condition,
                "passed": passed,
                "max_mass_balance_residual": mass_max,
                "max_bound_violation": bound_max,
                "max_integrality_error": integer_error,
                "flux_on_unselected_indicator_count": off_count,
                "max_flux_on_unselected_indicator": float(np.max(off_flux, initial=0)),
            }
        )
    return {
        "status": "passed" if all(r["passed"] for r in rows) else "failed",
        "feasibility_tolerance": feasibility_tolerance,
        "activity_tolerance": activity_tolerance,
        "conditions": rows,
        "claim_limit": "Numerical feasibility only; input biology unvalidated.",
    }


def validate_receipt(
    path: Path,
    *,
    context_sha256: str,
    conditions: Sequence[str],
    kind: str,
    required_mip_gap: float = 1e-4,
) -> dict[str, Any]:
    """Require context, exact identities, solver gap, primal audit and intact artifacts."""
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"receipt is not an object: {path}")
    expected = {
        "status": "completed",
        "scientific_success": True,
        "context_sha256": context_sha256,
        "schema_version": f"metabolic_{kind}_checkpoint.v1",
    }
    issues = [f"invalid_{key}" for key, wanted in expected.items() if value.get(key) != wanted]
    instrumentation = value.get("instrumentation", {})
    telemetry = instrumentation.get("telemetry", {})
    issues.extend(optimization_issues(telemetry))
    requested = telemetry.get("requested_mip_gap")
    if not _finite(requested) or requested > required_mip_gap + 1e-12:
        issues.append("requested_gap_is_looser_than_contract")
    if instrumentation.get("summary_error") is not None:
        issues.append("summary_error")
    validation = instrumentation.get("primal_validation", {})
    if validation.get("status") != "passed":
        issues.append("primal_validation_missing_or_failed")
    audit_conditions = [r.get("condition") for r in validation.get("conditions", [])]
    if audit_conditions != list(conditions) or not all(
        r.get("passed") is True for r in validation.get("conditions", [])
    ):
        issues.append("primal_audit_condition_mismatch")
    for row in validation.get("conditions", []):
        for key, maximum in (
            ("max_mass_balance_residual", 1e-6),
            ("max_bound_violation", 1e-6),
            ("max_integrality_error", 1e-5),
            ("flux_on_unselected_indicator_count", 0),
        ):
            observed = row.get(key)
            if not _finite(observed) or not 0 <= observed <= maximum:
                issues.append(f"primal_metric_failed_{key}")
    if kind == "independent":
        summaries = [value.get("solution", {})]
        if len(conditions) != 1 or value.get("condition") != conditions[0]:
            issues.append("condition_mismatch")
    elif kind == "joint":
        result = value.get("result", {})
        summaries = result.get("joint", [])
        if result.get("conditions") != list(conditions) or value.get("sample_count") != len(
            conditions
        ):
            issues.append("joint_condition_mismatch")
    else:
        raise ValueError(f"invalid receipt kind: {kind}")
    if [s.get("condition") for s in summaries] != list(conditions):
        issues.append("summary_condition_mismatch")
    for summary in summaries:
        if summary.get("status") != "optimal":
            issues.append("summary_not_optimal")
        if not _finite(summary.get("problem_objective_value")):
            issues.append("nonfinite_summary_objective")
        elif _finite(telemetry.get("objective_value")) and not math.isclose(
            summary["problem_objective_value"],
            telemetry["objective_value"],
            rel_tol=1e-8,
            abs_tol=1e-6,
        ):
            issues.append("summary_objective_mismatch")
        fluxes = summary.get("nonzero_fluxes", [])
        if any(not _finite(x[1]) for x in fluxes):
            issues.append("nonfinite_summary_flux")
        flux_ids = [x[0] for x in fluxes]
        if len(flux_ids) != len(set(flux_ids)):
            issues.append("duplicate_summary_reactions")
        if set(flux_ids) != set(summary.get("active_by_flux", [])):
            issues.append("summary_flux_set_mismatch")
        if set(flux_ids) - set(summary.get("active_by_indicator", [])):
            issues.append("flux_indicator_disagreement")
    if instrumentation.get("summaries") != summaries:
        issues.append("instrumentation_summary_mismatch")
    artifacts = instrumentation.get("artifacts", {})
    hashes = instrumentation.get("artifact_sha256", {})
    for label in ("solution", "mip_start", "gurobi_log"):
        raw = artifacts.get(label)
        artifact = Path(raw) if isinstance(raw, str) and raw else None
        if artifact is None or not artifact.is_file() or artifact.stat().st_size == 0:
            issues.append(f"missing_artifact_{label}")
        elif hashes.get(label) != sha256(artifact):
            issues.append(f"artifact_hash_mismatch_{label}")
    if issues:
        raise ValueError(f"receipt rejected {path}: {', '.join(issues)}")
    return value


def require_independent_receipts(context_path: Path) -> list[dict[str, Any]]:
    context = json.loads(context_path.read_text())
    digest = sha256(context_path)
    conditions = context["conditions"]
    if len(conditions) != len(set(conditions)):
        raise ValueError("duplicate context conditions")
    return [
        validate_receipt(
            context_path.parent / "independent" / f"{i:03d}_{c}.json",
            context_sha256=digest,
            conditions=[c],
            kind="independent",
        )
        for i, c in enumerate(conditions)
    ]
