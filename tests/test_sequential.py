from __future__ import annotations

from dataclasses import dataclass

from hgsoc_corneto.metabolic.sequential import (
    CandidateConstraint,
    apply_sequential_constraints,
    generate_meeson_candidates,
)


@dataclass
class FakeReaction:
    id: str
    bounds: tuple[float, float]
    gene_reaction_rule: str = ""

    @property
    def reversibility(self) -> bool:
        return self.bounds[0] < 0 < self.bounds[1]


class FakeReactionCollection(list[FakeReaction]):
    def get_by_id(self, reaction_id: str) -> FakeReaction:
        return next(reaction for reaction in self if reaction.id == reaction_id)


@dataclass
class FakeSolution:
    growth: float
    status: str = "optimal"

    @property
    def objective_value(self) -> float:
        return self.growth

    @property
    def fluxes(self) -> dict[str, float]:
        return {"BIOMASS": self.growth}


class ParallelFakeModel:
    def __init__(self, upper: float = 1000.0) -> None:
        self.reactions = FakeReactionCollection(
            [
                FakeReaction("A", (0.0, upper), "gene_a"),
                FakeReaction("B", (0.0, upper), "gene_b"),
            ]
        )

    def optimize(self) -> FakeSolution:
        total_routes = sum(reaction.bounds[1] for reaction in self.reactions)
        return FakeSolution(min(10.0, total_routes))


def candidate(reaction_id: str) -> CandidateConstraint:
    return CandidateConstraint(
        reaction_id=reaction_id,
        category="single_gene_forward",
        genes=(f"gene_{reaction_id.casefold()}",),
        expression_bound=2.0,
        proposed_lower=0.0,
        proposed_upper=2.0,
        reversible=False,
    )


def test_sequential_algorithm_is_order_dependent() -> None:
    first = apply_sequential_constraints(
        ParallelFakeModel(),
        [candidate("A"), candidate("B")],
        biomass_id="BIOMASS",
        growth_threshold=6.0,
        semantics="bounds_safe",
    )
    second = apply_sequential_constraints(
        ParallelFakeModel(),
        [candidate("B"), candidate("A")],
        biomass_id="BIOMASS",
        growth_threshold=6.0,
        semantics="bounds_safe",
    )

    assert first.retained_reactions == ["A"]
    assert first.reopened_reactions == ["B"]
    assert second.retained_reactions == ["B"]
    assert second.reopened_reactions == ["A"]
    assert first.final_growth == second.final_growth == 10.0


def test_published_reopen_can_expand_original_bound() -> None:
    model = ParallelFakeModel(upper=5.0)
    result = apply_sequential_constraints(
        model,
        [candidate("A"), candidate("B")],
        biomass_id="BIOMASS",
        growth_threshold=6.0,
        semantics="published",
    )

    reopened = result.decisions[1]
    assert reopened.previous_bounds == (0.0, 5.0)
    assert reopened.final_bounds == (0.0, 1000.0)


def test_bounds_safe_restores_original_bound() -> None:
    model = ParallelFakeModel(upper=5.0)
    result = apply_sequential_constraints(
        model,
        [candidate("A"), candidate("B")],
        biomass_id="BIOMASS",
        growth_threshold=6.0,
        semantics="bounds_safe",
    )

    reopened = result.decisions[1]
    assert reopened.final_bounds == (0.0, 5.0)


def test_candidate_generation_matches_upstream_group_order() -> None:
    model = type("Model", (), {})()
    model.reactions = FakeReactionCollection(
        [
            FakeReaction("ONE_F", (0.0, 1000.0), "g1"),
            FakeReaction("ONE_R", (-1000.0, 1000.0), "g2"),
            FakeReaction("OR_F", (0.0, 1000.0), "g3 or g4"),
            FakeReaction("AND_F", (0.0, 1000.0), "g5 and g6"),
            FakeReaction("MIXED", (0.0, 1000.0), "(g7 and g8) or g9"),
            FakeReaction("NONE", (0.0, 1000.0), ""),
        ]
    )

    candidates, audit = generate_meeson_candidates(
        model,
        {"g1": 1.0, "g2": 2.0, "g3": 3.0, "g4": 4.0, "g5": 5.0, "g6": 6.0},
    )

    assert [item.reaction_id for item in candidates] == ["ONE_F", "ONE_R", "OR_F", "AND_F"]
    assert [item.expression_bound for item in candidates] == [1.0, 2.0, 7.0, 5.0]
    assert audit["reaction_rule_partition"] == {
        "and": 1,
        "and_or": 1,
        "or": 1,
        "single_gene": 2,
        "no_gene": 1,
    }
    assert audit["skipped_reactions"]["and_or_not_implemented_upstream"] == ["MIXED"]
