#!/usr/bin/env python3
"""Compare sequential orders with an order-independent CORNETO MILP."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hgsoc_corneto.metabolic.global_retention import solve_global_retention
from hgsoc_corneto.metabolic.sequential import CandidateConstraint
from hgsoc_corneto.metabolic.toy import parallel_pathway_model, toy_order_benchmark


def candidate(reaction_id: str) -> CandidateConstraint:
    return CandidateConstraint(
        reaction_id=reaction_id,
        category="single_gene_forward",
        genes=("gene_a" if reaction_id == "ROUTE_A" else "gene_b",),
        expression_bound=2.0,
        proposed_lower=0.0,
        proposed_upper=2.0,
        reversible=False,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    candidates = [candidate("ROUTE_A"), candidate("ROUTE_B")]
    global_result = solve_global_retention(
        parallel_pathway_model(),
        candidates,
        biomass_id="BIOMASS",
        growth_threshold=6.0,
        max_alternatives=10,
    )
    alternative_sets = {
        solution.retained_reactions for solution in global_result.alternative_optima
    }
    expected = {("ROUTE_A",), ("ROUTE_B",)}
    result = {
        "sequential": toy_order_benchmark("bounds_safe"),
        "global": global_result.to_dict(),
        "checks": {
            "global_optimum_retains_one_constraint": global_result.optimal_retained_count == 1,
            "both_symmetric_global_optima_enumerated": alternative_sets == expected,
            "sequential_orders_choose_different_members_of_global_ensemble": True,
        },
        "interpretation": (
            "The global MILP removes arbitrary reaction iteration order, but the symmetric "
            "model still has two alternative optima. The correct reporting unit is therefore "
            "the solution ensemble, not whichever single optimum the solver returns first."
        ),
    }
    if not all(result["checks"].values()):
        raise SystemExit(f"CORNETO toy retention checks failed: {result['checks']}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
