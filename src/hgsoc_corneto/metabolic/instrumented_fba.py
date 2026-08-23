"""Time-bounded CORNETO sparse-FBA with auditable Gurobi telemetry."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from hgsoc_corneto.metabolic.joint_fba import (
    ConditionFluxSolution,
    _active_union,
    _corneto_components,
    _summaries,
)


@dataclass(frozen=True)
class SolverTelemetry:
    cvxpy_status: str
    gurobi_status_code: int | None
    gurobi_status_name: str | None
    has_incumbent: bool
    solution_count: int
    objective_value: float | None
    best_bound: float | None
    absolute_gap: float | None
    relative_gap: float | None
    node_count: float | None
    runtime_seconds: float | None
    work_units: float | None
    simplex_iterations: float | None
    barrier_iterations: int | None
    variable_count: int | None
    binary_variable_count: int | None
    integer_variable_count: int | None
    constraint_count: int | None
    nonzero_count: int | None
    requested_time_limit_seconds: int
    requested_mip_gap: float
    requested_threads: int


@dataclass(frozen=True)
class InstrumentedSparseFBAResult:
    status: str
    scientific_success: bool
    telemetry: SolverTelemetry
    summaries: tuple[ConditionFluxSolution, ...]
    active_union: tuple[str, ...]
    summary_error: str | None
    artifacts: Mapping[str, str]

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["active_union_size"] = len(self.active_union)
        return value


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _attribute(model: Any, name: str) -> Any:
    try:
        return getattr(model, name)
    except Exception:
        return None


def _status_name(code: int | None) -> str | None:
    if code is None:
        return None
    try:
        from gurobipy import GRB
    except ImportError:  # pragma: no cover - exercised on Roihu
        return None
    names = (
        "LOADED",
        "OPTIMAL",
        "INFEASIBLE",
        "INF_OR_UNBD",
        "UNBOUNDED",
        "CUTOFF",
        "ITERATION_LIMIT",
        "NODE_LIMIT",
        "TIME_LIMIT",
        "SOLUTION_LIMIT",
        "INTERRUPTED",
        "NUMERIC",
        "SUBOPTIMAL",
        "INPROGRESS",
        "USER_OBJ_LIMIT",
        "WORK_LIMIT",
        "MEM_LIMIT",
    )
    return next((name for name in names if getattr(GRB, name, None) == code), f"STATUS_{code}")


def _telemetry(
    solution: Any,
    *,
    max_seconds: int,
    mip_gap: float,
    threads: int,
) -> tuple[SolverTelemetry, Any]:
    stats = getattr(solution, "solver_stats", None)
    model = getattr(stats, "extra_stats", None)
    code_raw = _attribute(model, "Status")
    code = int(code_raw) if code_raw is not None else None
    count_raw = _attribute(model, "SolCount")
    solution_count = int(count_raw) if count_raw is not None else 0
    has_incumbent = solution_count > 0
    objective = _finite(_attribute(model, "ObjVal")) if has_incumbent else None
    bound = _finite(_attribute(model, "ObjBound"))
    relative_gap = _finite(_attribute(model, "MIPGap")) if has_incumbent else None
    absolute_gap = (
        abs(objective - bound) if objective is not None and bound is not None else None
    )
    telemetry = SolverTelemetry(
        cvxpy_status=str(getattr(solution, "status", "unknown")),
        gurobi_status_code=code,
        gurobi_status_name=_status_name(code),
        has_incumbent=has_incumbent,
        solution_count=solution_count,
        objective_value=objective,
        best_bound=bound,
        absolute_gap=absolute_gap,
        relative_gap=relative_gap,
        node_count=_finite(_attribute(model, "NodeCount")),
        runtime_seconds=_finite(_attribute(model, "Runtime")),
        work_units=_finite(_attribute(model, "Work")),
        simplex_iterations=_finite(_attribute(model, "IterCount")),
        barrier_iterations=(
            int(_attribute(model, "BarIterCount"))
            if _attribute(model, "BarIterCount") is not None
            else None
        ),
        variable_count=(
            int(_attribute(model, "NumVars")) if _attribute(model, "NumVars") is not None else None
        ),
        binary_variable_count=(
            int(_attribute(model, "NumBinVars"))
            if _attribute(model, "NumBinVars") is not None
            else None
        ),
        integer_variable_count=(
            int(_attribute(model, "NumIntVars"))
            if _attribute(model, "NumIntVars") is not None
            else None
        ),
        constraint_count=(
            int(_attribute(model, "NumConstrs"))
            if _attribute(model, "NumConstrs") is not None
            else None
        ),
        nonzero_count=(
            int(_attribute(model, "NumNZs")) if _attribute(model, "NumNZs") is not None else None
        ),
        requested_time_limit_seconds=max_seconds,
        requested_mip_gap=mip_gap,
        requested_threads=threads,
    )
    return telemetry, model


def _write_incumbent_artifacts(model: Any, prefix: Path, has_incumbent: bool) -> dict[str, str]:
    artifacts: dict[str, str] = {}
    if model is None or not has_incumbent:
        return artifacts
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for label, suffix in (("solution", ".sol"), ("mip_start", ".mst")):
        path = prefix.with_suffix(suffix)
        try:
            model.write(str(path))
        except Exception as error:  # pragma: no cover - solver/filesystem specific
            artifacts[f"{label}_error"] = f"{type(error).__name__}: {error}"
        else:
            artifacts[label] = str(path)
    return artifacts


def _finish(
    problem: Any,
    solution: Any,
    reaction_ids: Sequence[str],
    conditions: Sequence[str],
    *,
    max_seconds: int,
    mip_gap: float,
    threads: int,
    active_tolerance: float,
    artifact_prefix: Path,
    log_file: Path,
) -> InstrumentedSparseFBAResult:
    telemetry, model = _telemetry(
        solution, max_seconds=max_seconds, mip_gap=mip_gap, threads=threads
    )
    summaries: tuple[ConditionFluxSolution, ...] = ()
    summary_error = None
    if telemetry.has_incumbent:
        try:
            summaries = _summaries(
                problem,
                solution,
                reaction_ids,
                conditions,
                active_tolerance=active_tolerance,
            )
        except Exception as error:  # preserve telemetry even if expression extraction fails
            summary_error = f"{type(error).__name__}: {error}"
    normalized = telemetry.cvxpy_status.casefold()
    scientific_success = normalized in {"optimal", "optimal_inaccurate"} and bool(summaries)
    if scientific_success:
        status = "completed"
    elif telemetry.has_incumbent:
        status = "partial_incumbent"
    elif telemetry.gurobi_status_name == "TIME_LIMIT":
        status = "time_limit_no_incumbent"
    else:
        status = "failed_no_incumbent"
    artifacts = {"gurobi_log": str(log_file)}
    artifacts.update(_write_incumbent_artifacts(model, artifact_prefix, telemetry.has_incumbent))
    return InstrumentedSparseFBAResult(
        status=status,
        scientific_success=scientific_success,
        telemetry=telemetry,
        summaries=summaries,
        active_union=_active_union(summaries),
        summary_error=summary_error,
        artifacts=artifacts,
    )


def solve_independent_instrumented(
    model: Any,
    *,
    condition: str,
    objective: Mapping[str, float],
    reaction_bounds: Mapping[str, tuple[float | None, float | None]],
    independent_lambda: float,
    max_seconds: int,
    mip_gap: float,
    threads: int,
    artifact_prefix: Path,
    log_file: Path,
    solver: str = "gurobi",
    active_tolerance: float = 1e-7,
) -> InstrumentedSparseFBAResult:
    """Solve one condition while preserving a time-limited incumbent and gap."""

    MultiSampleFBA, graph, reaction_ids = _corneto_components(model)
    problem = MultiSampleFBA(lambda_reg=independent_lambda).build(
        graph,
        objectives=dict(objective),
        reaction_bounds=dict(reaction_bounds),
    )
    log_file.parent.mkdir(parents=True, exist_ok=True)
    solution = problem.solve(
        solver=solver,
        max_seconds=max_seconds,
        verbosity=1,
        MIPGap=mip_gap,
        Threads=threads,
        Seed=0,
        LogFile=str(log_file),
    )
    return _finish(
        problem,
        solution,
        reaction_ids,
        (condition,),
        max_seconds=max_seconds,
        mip_gap=mip_gap,
        threads=threads,
        active_tolerance=active_tolerance,
        artifact_prefix=artifact_prefix,
        log_file=log_file,
    )


def solve_joint_instrumented(
    model: Any,
    *,
    objectives: Mapping[str, Mapping[str, float]],
    reaction_bounds: Mapping[str, Mapping[str, tuple[float | None, float | None]]],
    joint_lambda: float,
    max_seconds: int,
    mip_gap: float,
    threads: int,
    artifact_prefix: Path,
    log_file: Path,
    solver: str = "gurobi",
    active_tolerance: float = 1e-7,
) -> InstrumentedSparseFBAResult:
    """Solve a cohort jointly while preserving a time-limited incumbent and gap."""

    conditions = tuple(objectives)
    MultiSampleFBA, graph, reaction_ids = _corneto_components(model)
    problem = MultiSampleFBA(lambda_reg=joint_lambda).build_many(
        graph,
        objectives={condition: dict(objectives[condition]) for condition in conditions},
        reaction_bounds={condition: dict(reaction_bounds[condition]) for condition in conditions},
    )
    log_file.parent.mkdir(parents=True, exist_ok=True)
    solution = problem.solve(
        solver=solver,
        max_seconds=max_seconds,
        verbosity=1,
        MIPGap=mip_gap,
        Threads=threads,
        Seed=0,
        LogFile=str(log_file),
    )
    return _finish(
        problem,
        solution,
        reaction_ids,
        conditions,
        max_seconds=max_seconds,
        mip_gap=mip_gap,
        threads=threads,
        active_tolerance=active_tolerance,
        artifact_prefix=artifact_prefix,
        log_file=log_file,
    )
