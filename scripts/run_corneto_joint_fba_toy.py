#!/usr/bin/env python3
"""Demonstrate CORNETO joint structured sparsity on two toy samples."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hgsoc_corneto.metabolic.joint_fba import compare_independent_and_joint_sparse_fba
from hgsoc_corneto.metabolic.toy import parallel_pathway_model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    objectives = {
        "sample_prefers_A": {"ROUTE_B": 1.0},
        "sample_prefers_B": {"ROUTE_A": 1.0},
    }
    reaction_bounds = {
        condition: {"BIOMASS": (6.0, 6.0)} for condition in objectives
    }
    comparison = compare_independent_and_joint_sparse_fba(
        parallel_pathway_model(),
        objectives=objectives,
        reaction_bounds=reaction_bounds,
        independent_lambda=0.1,
        joint_lambda=10.0,
    )
    independent_routes = {
        solution.condition: tuple(
            reaction
            for reaction in solution.active_by_flux
            if reaction in {"ROUTE_A", "ROUTE_B"}
        )
        for solution in comparison.independent
    }
    joint_routes = {
        solution.condition: tuple(
            reaction
            for reaction in solution.active_by_flux
            if reaction in {"ROUTE_A", "ROUTE_B"}
        )
        for solution in comparison.joint
    }
    joint_route_union = {reaction for routes in joint_routes.values() for reaction in routes}
    checks = {
        "independent_samples_choose_opposite_routes": independent_routes
        == {
            "sample_prefers_A": ("ROUTE_A",),
            "sample_prefers_B": ("ROUTE_B",),
        },
        "joint_samples_share_one_route": len(joint_route_union) == 1
        and len(set(joint_routes.values())) == 1,
        "joint_active_union_is_smaller": len(comparison.joint_active_union)
        < len(comparison.independent_active_union),
    }
    result = {
        "design": {
            "model": "two interchangeable routes with shared uptake",
            "biomass_fixed_per_sample": 6.0,
            "sample_specific_linear_cost": 1.0,
            "independent_lambda": 0.1,
            "joint_lambda": 10.0,
        },
        "comparison": comparison.to_dict(),
        "independent_route_choices": independent_routes,
        "joint_route_choices": joint_routes,
        "checks": checks,
        "interpretation": (
            "Independent sparse FBA follows each sample's weak route preference, whereas "
            "CORNETO's joint logical-OR penalty selects a shared route and reduces the "
            "cross-sample active-reaction union. This verifies the multi-sample mechanism; "
            "it is not yet evidence from OCM data."
        ),
    }
    if not all(checks.values()):
        raise SystemExit(f"Joint sparse-FBA toy checks failed: {checks}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
