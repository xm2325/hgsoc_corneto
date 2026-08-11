from __future__ import annotations

import importlib.util

import pytest

from hgsoc_corneto.metabolic.global_retention import solve_global_retention
from hgsoc_corneto.metabolic.sequential import CandidateConstraint
from hgsoc_corneto.metabolic.toy import parallel_pathway_model

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("corneto") is None or importlib.util.find_spec("cobra") is None,
    reason="CORNETO metabolic tests run in the pinned Roihu environment",
)


def test_global_retention_enumerates_symmetric_optima() -> None:
    candidates = [
        CandidateConstraint(
            reaction_id=reaction_id,
            category="single_gene_forward",
            genes=(gene,),
            expression_bound=2.0,
            proposed_lower=0.0,
            proposed_upper=2.0,
            reversible=False,
        )
        for reaction_id, gene in (("ROUTE_A", "gene_a"), ("ROUTE_B", "gene_b"))
    ]

    result = solve_global_retention(
        parallel_pathway_model(),
        candidates,
        biomass_id="BIOMASS",
        growth_threshold=6.0,
        max_alternatives=10,
    )

    assert result.optimal_retained_count == 1
    assert {solution.retained_reactions for solution in result.alternative_optima} == {
        ("ROUTE_A",),
        ("ROUTE_B",),
    }
