"""Auditable reconstruction of the sequential Meeson constraint algorithm."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from typing import Any, Literal

Semantics = Literal["published", "bounds_safe"]


@dataclass(frozen=True)
class CandidateConstraint:
    reaction_id: str
    category: str
    genes: tuple[str, ...]
    expression_bound: float
    proposed_lower: float
    proposed_upper: float
    reversible: bool

    @property
    def proposed_bounds(self) -> tuple[float, float]:
        return (self.proposed_lower, self.proposed_upper)


@dataclass(frozen=True)
class ConstraintDecision:
    order_index: int
    reaction_id: str
    category: str
    genes: tuple[str, ...]
    previous_bounds: tuple[float, float]
    proposed_bounds: tuple[float, float]
    trial_bounds: tuple[float, float]
    trial_status: str
    trial_growth: float | None
    retained: bool
    final_bounds: tuple[float, float]
    post_decision_status: str
    post_decision_growth: float | None
    reason: str


@dataclass
class SequentialResult:
    semantics: Semantics
    growth_threshold: float
    initial_status: str
    initial_growth: float | None
    final_status: str
    final_growth: float | None
    decisions: list[ConstraintDecision]

    @property
    def retained_reactions(self) -> list[str]:
        return [decision.reaction_id for decision in self.decisions if decision.retained]

    @property
    def reopened_reactions(self) -> list[str]:
        return [decision.reaction_id for decision in self.decisions if not decision.retained]

    def to_dict(self) -> dict[str, Any]:
        return {
            "semantics": self.semantics,
            "growth_threshold": self.growth_threshold,
            "initial_status": self.initial_status,
            "initial_growth": self.initial_growth,
            "final_status": self.final_status,
            "final_growth": self.final_growth,
            "retained_reactions": self.retained_reactions,
            "reopened_reactions": self.reopened_reactions,
            "decisions": [asdict(decision) for decision in self.decisions],
        }


def _solution_status(solution: Any) -> str:
    return str(getattr(solution, "status", "unknown"))


def _solution_growth(solution: Any, biomass_id: str) -> float | None:
    if _solution_status(solution).casefold() != "optimal":
        return None
    try:
        value = float(solution.fluxes[biomass_id])
    except (AttributeError, KeyError, TypeError, ValueError):
        value = float(getattr(solution, "objective_value", math.nan))
    return value if math.isfinite(value) else None


def _trial_bounds(
    previous: tuple[float, float], proposed: tuple[float, float], semantics: Semantics
) -> tuple[float, float]:
    if semantics == "published":
        lower, upper = proposed
    elif semantics == "bounds_safe":
        lower = max(previous[0], proposed[0])
        upper = min(previous[1], proposed[1])
    else:  # pragma: no cover - protected by the public type and explicit validation
        raise ValueError(f"Unknown semantics: {semantics}")
    if not (math.isfinite(lower) and math.isfinite(upper)):
        raise ValueError(f"Non-finite candidate bounds: {(lower, upper)}")
    if lower > upper:
        raise ValueError(f"Invalid candidate bounds: {(lower, upper)}")
    return (lower, upper)


def _reopen_bounds(
    previous: tuple[float, float], candidate: CandidateConstraint, semantics: Semantics
) -> tuple[float, float]:
    if semantics == "bounds_safe":
        return previous
    return (-1000.0, 1000.0) if candidate.reversible else (0.0, 1000.0)


def apply_sequential_constraints(
    model: Any,
    candidates: Iterable[CandidateConstraint],
    *,
    biomass_id: str,
    growth_threshold: float,
    semantics: Semantics = "published",
) -> SequentialResult:
    """Mutate ``model`` by testing candidate bounds in the supplied order.

    ``published`` reproduces the public implementation's hard-coded reopening
    bounds. ``bounds_safe`` restores each reaction's pre-trial bounds and never
    expands the feasible interval when expression exceeds the original bound.
    Both retain a constraint only when optimized growth is strictly greater than
    the experimental threshold, matching the public ``<=`` reopening test.
    """

    if semantics not in {"published", "bounds_safe"}:
        raise ValueError(f"Unknown semantics: {semantics}")
    if not math.isfinite(growth_threshold) or growth_threshold < 0:
        raise ValueError("growth_threshold must be finite and non-negative")

    initial_solution = model.optimize()
    initial_status = _solution_status(initial_solution)
    initial_growth = _solution_growth(initial_solution, biomass_id)
    decisions: list[ConstraintDecision] = []

    for order_index, candidate in enumerate(candidates):
        reaction = model.reactions.get_by_id(candidate.reaction_id)
        previous = tuple(float(value) for value in reaction.bounds)
        trial = _trial_bounds(previous, candidate.proposed_bounds, semantics)
        reaction.bounds = trial
        trial_solution = model.optimize()
        trial_status = _solution_status(trial_solution)
        trial_growth = _solution_growth(trial_solution, biomass_id)

        retained = trial_growth is not None and trial_growth > growth_threshold
        if retained:
            final_bounds = trial
            post_solution = trial_solution
            reason = "growth_strictly_above_threshold"
        else:
            final_bounds = _reopen_bounds(previous, candidate, semantics)
            reaction.bounds = final_bounds
            post_solution = model.optimize()
            reason = (
                "trial_not_optimal"
                if trial_growth is None
                else "growth_at_or_below_threshold"
            )

        decisions.append(
            ConstraintDecision(
                order_index=order_index,
                reaction_id=candidate.reaction_id,
                category=candidate.category,
                genes=candidate.genes,
                previous_bounds=previous,
                proposed_bounds=candidate.proposed_bounds,
                trial_bounds=trial,
                trial_status=trial_status,
                trial_growth=trial_growth,
                retained=retained,
                final_bounds=final_bounds,
                post_decision_status=_solution_status(post_solution),
                post_decision_growth=_solution_growth(post_solution, biomass_id),
                reason=reason,
            )
        )

    final_solution = model.optimize()
    return SequentialResult(
        semantics=semantics,
        growth_threshold=growth_threshold,
        initial_status=initial_status,
        initial_growth=initial_growth,
        final_status=_solution_status(final_solution),
        final_growth=_solution_growth(final_solution, biomass_id),
        decisions=decisions,
    )


def _published_gene_tokens(rule: str) -> tuple[str, ...]:
    tokens = []
    for token in rule.split():
        if token not in {"and", "or", "(", ")"}:
            gene = token.strip("()")
            if gene:
                tokens.append(gene)
    return tuple(tokens)


def _classify_rule(rule: str) -> str:
    has_and = "and" in rule
    has_or = "or" in rule
    if has_and and has_or:
        return "and_or"
    if has_and:
        return "and"
    if has_or:
        return "or"
    if not rule:
        return "no_gene"
    return "single_gene"


def generate_meeson_candidates(
    model: Any,
    expression: Mapping[str, float],
    *,
    media_reactions: Iterable[str] = (),
    essential_reactions: Iterable[str] = (),
    essential_genes: Iterable[str] = (),
) -> tuple[list[CandidateConstraint], dict[str, Any]]:
    """Generate candidates in the exact category order used upstream.

    Mixed AND/OR reactions are reported but omitted because the public function
    classifies them and never supplies an integration loop for them.
    """

    normalized_expression = {str(gene): float(value) for gene, value in expression.items()}
    invalid = {
        gene: value
        for gene, value in normalized_expression.items()
        if not math.isfinite(value) or value < 0
    }
    if invalid:
        raise ValueError(f"Expression bounds must be finite and non-negative: {invalid}")

    media = set(media_reactions)
    essential_rxns = set(essential_reactions)
    essential_gene_ids = set(essential_genes)
    grouped: dict[str, list[Any]] = {
        "and": [],
        "and_or": [],
        "or": [],
        "single_gene": [],
        "no_gene": [],
    }
    for reaction in model.reactions:
        grouped[_classify_rule(reaction.gene_reaction_rule)].append(reaction)

    ordered_groups = (
        ("single_gene_forward", grouped["single_gene"], False),
        ("single_gene_reversible", grouped["single_gene"], True),
        ("or_forward", grouped["or"], False),
        ("or_reversible", grouped["or"], True),
        ("and_forward", grouped["and"], False),
        ("and_reversible", grouped["and"], True),
    )
    candidates: list[CandidateConstraint] = []
    skipped: dict[str, list[str]] = {
        "media": [],
        "essential_reaction": [],
        "essential_single_gene": [],
        "single_gene_absent_from_expression": [],
        "and_or_not_implemented_upstream": [reaction.id for reaction in grouped["and_or"]],
        "no_gene": [reaction.id for reaction in grouped["no_gene"]],
    }

    for category, reactions, reversible in ordered_groups:
        for reaction in reactions:
            if bool(reaction.reversibility) is not reversible:
                continue
            if reaction.id in media:
                skipped["media"].append(reaction.id)
                continue
            if reaction.id in essential_rxns:
                skipped["essential_reaction"].append(reaction.id)
                continue

            rule = reaction.gene_reaction_rule
            genes = _published_gene_tokens(rule)
            if category.startswith("single_gene"):
                if rule in essential_gene_ids:
                    skipped["essential_single_gene"].append(reaction.id)
                    continue
                if rule not in normalized_expression:
                    skipped["single_gene_absent_from_expression"].append(reaction.id)
                    continue
                expression_bound = normalized_expression[rule]
            else:
                present = [
                    normalized_expression[gene]
                    for gene in genes
                    if gene in normalized_expression
                ]
                expression_bound = (
                    sum(present)
                    if category.startswith("or")
                    else min(present, default=0.0)
                )

            lower = -expression_bound if reversible else 0.0
            candidates.append(
                CandidateConstraint(
                    reaction_id=reaction.id,
                    category=category,
                    genes=genes,
                    expression_bound=expression_bound,
                    proposed_lower=lower,
                    proposed_upper=expression_bound,
                    reversible=reversible,
                )
            )

    audit = {
        "reaction_rule_partition": {key: len(value) for key, value in grouped.items()},
        "candidate_count": len(candidates),
        "candidate_counts_by_category": {
            category: sum(candidate.category == category for candidate in candidates)
            for category, _, _ in ordered_groups
        },
        "skipped_counts": {key: len(value) for key, value in skipped.items()},
        "skipped_reactions": skipped,
    }
    return candidates, audit
