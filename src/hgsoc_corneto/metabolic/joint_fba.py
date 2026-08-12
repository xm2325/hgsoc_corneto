"""Auditable independent-versus-joint CORNETO sparse-FBA comparison."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class ConditionFluxSolution:
    condition: str
    status: str
    problem_objective_value: float | None
    active_by_flux: tuple[str, ...]
    active_by_indicator: tuple[str, ...]
    nonzero_fluxes: tuple[tuple[str, float], ...]


@dataclass(frozen=True)
class SparseFBAComparison:
    solver: str
    independent_lambda: float
    joint_lambda: float
    active_tolerance: float
    conditions: tuple[str, ...]
    independent: tuple[ConditionFluxSolution, ...]
    joint: tuple[ConditionFluxSolution, ...]
    independent_active_union: tuple[str, ...]
    joint_active_union: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["independent_active_union_size"] = len(self.independent_active_union)
        result["joint_active_union_size"] = len(self.joint_active_union)
        return result


@dataclass(frozen=True)
class JointSparseFBAResult:
    """Joint-only sparse-FBA result for a pre-defined condition collection."""

    solver: str
    joint_lambda: float
    active_tolerance: float
    conditions: tuple[str, ...]
    joint: tuple[ConditionFluxSolution, ...]
    joint_active_union: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["joint_active_union_size"] = len(self.joint_active_union)
        return result


def _corneto_components(model: Any) -> tuple[Any, Any, list[str]]:
    try:
        from corneto.io import cobra_model_to_graph
        from corneto.methods.fba import MultiSampleFBA
    except ImportError as error:  # pragma: no cover - exercised on Roihu
        raise RuntimeError("Install pinned CORNETO and the 'metabolic' extra") from error

    graph = cobra_model_to_graph(model)
    reaction_ids = [str(value) for value in graph.get_attr_from_edges("id")]
    if len(reaction_ids) != len(model.reactions):
        raise ValueError("COBRA-to-CORNETO reaction ordering is inconsistent")
    return MultiSampleFBA, graph, reaction_ids


def _objective_value(solution: Any) -> float | None:
    raw = getattr(solution, "value", None)
    if raw is None:
        return None
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def _summaries(
    problem: Any,
    solution: Any,
    reaction_ids: Sequence[str],
    conditions: Sequence[str],
    *,
    active_tolerance: float,
) -> tuple[ConditionFluxSolution, ...]:
    import numpy as np

    flows = np.asarray(problem.expr.flow.value, dtype=float)
    indicators = np.asarray(problem.expr.edge_has_flux.value, dtype=float)
    if len(conditions) == 1:
        flows = flows.reshape(-1, 1)
        indicators = indicators.reshape(-1, 1)
    expected = (len(reaction_ids), len(conditions))
    if flows.shape != expected or indicators.shape != expected:
        raise ValueError(
            f"Unexpected CORNETO solution shapes: flow={flows.shape}, "
            f"indicator={indicators.shape}, expected={expected}"
        )

    status = str(getattr(solution, "status", "unknown"))
    objective = _objective_value(solution)
    summaries = []
    for column, condition in enumerate(conditions):
        nonzero = tuple(
            (reaction_id, float(flows[row, column]))
            for row, reaction_id in enumerate(reaction_ids)
            if abs(flows[row, column]) > active_tolerance
        )
        active_by_flux = tuple(reaction_id for reaction_id, _ in nonzero)
        active_by_indicator = tuple(
            reaction_id
            for row, reaction_id in enumerate(reaction_ids)
            if indicators[row, column] >= 0.5
        )
        summaries.append(
            ConditionFluxSolution(
                condition=condition,
                status=status,
                problem_objective_value=objective,
                active_by_flux=active_by_flux,
                active_by_indicator=active_by_indicator,
                nonzero_fluxes=nonzero,
            )
        )
    return tuple(summaries)


def _active_union(solutions: Sequence[ConditionFluxSolution]) -> tuple[str, ...]:
    active = {
        reaction for solution in solutions for reaction in solution.active_by_flux
    }
    return tuple(sorted(active))


def solve_joint_sparse_fba(
    model: Any,
    *,
    objectives: Mapping[str, Mapping[str, float]],
    reaction_bounds: Mapping[
        str,
        Mapping[str, tuple[float | None, float | None]],
    ],
    joint_lambda: float,
    solver: str = "highs",
    active_tolerance: float = 1e-7,
) -> JointSparseFBAResult:
    """Solve one joint union-sparse FBA without repeating independent MILPs."""

    conditions = tuple(objectives)
    if not conditions:
        raise ValueError("At least one condition is required")
    if set(conditions) != set(reaction_bounds):
        raise ValueError("objectives and reaction_bounds must name the same conditions")
    if joint_lambda < 0:
        raise ValueError("Regularization weight must be non-negative")
    if active_tolerance <= 0:
        raise ValueError("active_tolerance must be positive")

    MultiSampleFBA, graph, reaction_ids = _corneto_components(model)
    problem = MultiSampleFBA(lambda_reg=joint_lambda).build_many(
        graph,
        objectives={condition: dict(objectives[condition]) for condition in conditions},
        reaction_bounds={
            condition: dict(reaction_bounds[condition]) for condition in conditions
        },
    )
    solution = problem.solve(solver=solver)
    summaries = _summaries(
        problem,
        solution,
        reaction_ids,
        conditions,
        active_tolerance=active_tolerance,
    )
    if any(
        summary.status.casefold() not in {"optimal", "optimal_inaccurate"}
        for summary in summaries
    ):
        raise RuntimeError(f"Joint conditions failed: {summaries[0].status}")
    return JointSparseFBAResult(
        solver=solver,
        joint_lambda=joint_lambda,
        active_tolerance=active_tolerance,
        conditions=conditions,
        joint=summaries,
        joint_active_union=_active_union(summaries),
    )


def compare_independent_and_joint_sparse_fba(
    model: Any,
    *,
    objectives: Mapping[str, Mapping[str, float]],
    reaction_bounds: Mapping[
        str,
        Mapping[str, tuple[float | None, float | None]],
    ],
    independent_lambda: float,
    joint_lambda: float,
    solver: str = "highs",
    active_tolerance: float = 1e-7,
) -> SparseFBAComparison:
    """Solve the same conditions independently and with joint union sparsity.

    The independent problems penalize each condition's active reactions. The
    joint problem uses CORNETO's structured regularizer over the logical OR of
    reaction-activity indicators across conditions.
    """

    conditions = tuple(objectives)
    if not conditions:
        raise ValueError("At least one condition is required")
    if set(conditions) != set(reaction_bounds):
        raise ValueError("objectives and reaction_bounds must name the same conditions")
    if independent_lambda < 0 or joint_lambda < 0:
        raise ValueError("Regularization weights must be non-negative")
    if active_tolerance <= 0:
        raise ValueError("active_tolerance must be positive")

    MultiSampleFBA, graph, reaction_ids = _corneto_components(model)
    independent = []
    for condition in conditions:
        problem = MultiSampleFBA(lambda_reg=independent_lambda).build(
            graph,
            objectives=dict(objectives[condition]),
            reaction_bounds=dict(reaction_bounds[condition]),
        )
        solution = problem.solve(solver=solver)
        summary = _summaries(
            problem,
            solution,
            reaction_ids,
            (condition,),
            active_tolerance=active_tolerance,
        )[0]
        if summary.status.casefold() not in {"optimal", "optimal_inaccurate"}:
            raise RuntimeError(f"Independent condition {condition!r} failed: {summary.status}")
        independent.append(summary)

    joint_problem = MultiSampleFBA(lambda_reg=joint_lambda).build_many(
        graph,
        objectives={condition: dict(objectives[condition]) for condition in conditions},
        reaction_bounds={
            condition: dict(reaction_bounds[condition]) for condition in conditions
        },
    )
    joint_solution = joint_problem.solve(solver=solver)
    joint = _summaries(
        joint_problem,
        joint_solution,
        reaction_ids,
        conditions,
        active_tolerance=active_tolerance,
    )
    if any(
        summary.status.casefold() not in {"optimal", "optimal_inaccurate"}
        for summary in joint
    ):
        raise RuntimeError(f"Joint conditions failed: {joint[0].status}")

    independent_tuple = tuple(independent)
    return SparseFBAComparison(
        solver=solver,
        independent_lambda=independent_lambda,
        joint_lambda=joint_lambda,
        active_tolerance=active_tolerance,
        conditions=conditions,
        independent=independent_tuple,
        joint=joint,
        independent_active_union=_active_union(independent_tuple),
        joint_active_union=_active_union(joint),
    )
