"""Order-independent CORNETO formulation of expression-bound retention."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from hgsoc_corneto.metabolic.sequential import CandidateConstraint


@dataclass(frozen=True)
class RetentionSolution:
    status: str
    objective_value: float | None
    growth: float | None
    retained_reactions: tuple[str, ...]
    reopened_reactions: tuple[str, ...]
    retained_vector: tuple[int, ...]


@dataclass(frozen=True)
class GlobalRetentionResult:
    formulation: str
    solver: str
    biomass_id: str
    growth_threshold: float
    strict_margin: float
    candidate_reactions: tuple[str, ...]
    primary: RetentionSolution
    optimal_retained_count: int
    alternative_optima: tuple[RetentionSolution, ...]
    alternatives_truncated: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _corneto_components(model: Any) -> tuple[Any, Any, list[str], Any, Any]:
    try:
        import corneto as cn
        import numpy as np
        from corneto.io import cobra_model_to_graph
    except ImportError as error:  # pragma: no cover - exercised on Roihu
        raise RuntimeError("Install CORNETO and the 'metabolic' extra") from error

    graph = cobra_model_to_graph(model)
    reaction_ids = list(graph.get_attr_from_edges("id"))
    lower = np.asarray([reaction.lower_bound for reaction in model.reactions], dtype=float)
    upper = np.asarray([reaction.upper_bound for reaction in model.reactions], dtype=float)
    if len(reaction_ids) != len(lower):
        raise ValueError("COBRA-to-CORNETO reaction ordering is inconsistent")
    return cn, graph, reaction_ids, lower, upper


def _build_problem(
    model: Any,
    candidates: list[CandidateConstraint],
    *,
    biomass_id: str,
    growth_threshold: float,
    strict_margin: float,
) -> tuple[Any, Any, Any, int, list[str]]:
    cn, graph, reaction_ids, lower, upper = _corneto_components(model)
    reaction_index = {reaction_id: index for index, reaction_id in enumerate(reaction_ids)}
    if biomass_id not in reaction_index:
        raise KeyError(f"Unknown biomass reaction: {biomass_id}")
    missing = [
        candidate.reaction_id
        for candidate in candidates
        if candidate.reaction_id not in reaction_index
    ]
    if missing:
        raise KeyError(f"Unknown candidate reactions: {missing}")
    if len({candidate.reaction_id for candidate in candidates}) != len(candidates):
        raise ValueError("Each reaction may appear at most once in a retention problem")

    problem = cn.opt.Flow(graph, lb=lower, ub=upper, values=True)
    flow = problem.expr.flow
    retained = cn.opt.Variable(
        "constraint_retained",
        (len(candidates),),
        vartype=cn.VarType.BINARY,
    )
    biomass_index = reaction_index[biomass_id]
    problem += flow[biomass_index] >= growth_threshold + strict_margin

    for candidate_index, candidate in enumerate(candidates):
        index = reaction_index[candidate.reaction_id]
        original_lower = lower[index]
        original_upper = upper[index]
        candidate_lower = max(original_lower, candidate.proposed_lower)
        candidate_upper = min(original_upper, candidate.proposed_upper)
        if candidate_lower > candidate_upper:
            raise ValueError(
                f"Candidate {candidate.reaction_id} has an empty bounds intersection"
            )
        # z=1 applies the expression-derived interval; z=0 restores the
        # original model interval. All candidates are optimized simultaneously.
        problem += flow[index] >= original_lower + (
            candidate_lower - original_lower
        ) * retained[candidate_index]
        problem += flow[index] <= original_upper - (
            original_upper - candidate_upper
        ) * retained[candidate_index]

    problem.add_objective(-retained.sum(), name="maximize_retained_constraints")
    return problem, retained, flow, biomass_index, reaction_ids


def _extract_solution(
    backend_solution: Any,
    retained: Any,
    flow: Any,
    biomass_index: int,
    candidates: list[CandidateConstraint],
) -> RetentionSolution:
    import numpy as np

    status = str(getattr(backend_solution, "status", "unknown"))
    raw_objective = getattr(backend_solution, "value", None)
    objective = float(raw_objective) if raw_objective is not None else None
    if retained.value is None or flow.value is None:
        vector = tuple(0 for _ in candidates)
        growth = None
    else:
        vector = tuple(int(value >= 0.5) for value in np.asarray(retained.value).reshape(-1))
        growth = float(np.asarray(flow.value).reshape(-1)[biomass_index])
    kept = tuple(
        candidate.reaction_id
        for candidate, value in zip(candidates, vector, strict=True)
        if value
    )
    reopened = tuple(
        candidate.reaction_id
        for candidate, value in zip(candidates, vector, strict=True)
        if not value
    )
    return RetentionSolution(
        status=status,
        objective_value=objective,
        growth=growth,
        retained_reactions=kept,
        reopened_reactions=reopened,
        retained_vector=vector,
    )


def solve_global_retention(
    model: Any,
    candidates: list[CandidateConstraint],
    *,
    biomass_id: str,
    growth_threshold: float,
    strict_margin: float = 1e-6,
    solver: str = "highs",
    enumerate_alternatives: bool = True,
    max_alternatives: int = 100,
) -> GlobalRetentionResult:
    """Maximize retained expression constraints in a single global MILP.

    The strict inequality in the sequential implementation is represented as
    ``growth >= threshold + strict_margin``. Alternative optimal retained sets
    are enumerated with no-good cuts, up to ``max_alternatives``.
    """

    if not candidates:
        raise ValueError("At least one candidate constraint is required")
    if strict_margin <= 0:
        raise ValueError("strict_margin must be positive")
    if max_alternatives < 1:
        raise ValueError("max_alternatives must be at least one")

    problem, retained, flow, biomass_index, _ = _build_problem(
        model,
        candidates,
        biomass_id=biomass_id,
        growth_threshold=growth_threshold,
        strict_margin=strict_margin,
    )
    backend_solution = problem.solve(solver=solver)
    primary = _extract_solution(backend_solution, retained, flow, biomass_index, candidates)
    if primary.status.casefold() not in {"optimal", "optimal_inaccurate"}:
        raise RuntimeError(f"Global retention problem did not solve: {primary.status}")
    optimum_count = sum(primary.retained_vector)

    alternatives: list[RetentionSolution] = [primary]
    truncated = False
    if enumerate_alternatives:
        problem += retained.sum() == optimum_count
        while len(alternatives) < max_alternatives:
            vector = alternatives[-1].retained_vector
            ones = [index for index, value in enumerate(vector) if value]
            zeros = [index for index, value in enumerate(vector) if not value]
            difference = 0
            if ones:
                difference += (1 - retained[ones]).sum()
            if zeros:
                difference += retained[zeros].sum()
            problem += difference >= 1
            backend_solution = problem.solve(solver=solver)
            candidate_solution = _extract_solution(
                backend_solution, retained, flow, biomass_index, candidates
            )
            if candidate_solution.status.casefold() not in {"optimal", "optimal_inaccurate"}:
                break
            alternatives.append(candidate_solution)
        else:
            truncated = True

    return GlobalRetentionResult(
        formulation=(
            "CORNETO stoichiometric flow MILP; maximize simultaneous retention of "
            "expression-derived reaction bounds subject to experimental growth"
        ),
        solver=solver,
        biomass_id=biomass_id,
        growth_threshold=growth_threshold,
        strict_margin=strict_margin,
        candidate_reactions=tuple(candidate.reaction_id for candidate in candidates),
        primary=primary,
        optimal_retained_count=optimum_count,
        alternative_optima=tuple(alternatives),
        alternatives_truncated=truncated,
    )
