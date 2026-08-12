from __future__ import annotations

import importlib.util

import pytest

from hgsoc_corneto.metabolic.joint_fba import (
    compare_independent_and_joint_sparse_fba,
    solve_joint_sparse_fba,
)
from hgsoc_corneto.metabolic.toy import parallel_pathway_model

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("corneto") is None or importlib.util.find_spec("cobra") is None,
    reason="CORNETO metabolic tests run in the pinned Roihu environment",
)


def test_joint_sparse_fba_reduces_cross_sample_union() -> None:
    objectives = {
        "sample_prefers_A": {"ROUTE_B": 1.0},
        "sample_prefers_B": {"ROUTE_A": 1.0},
    }
    bounds = {condition: {"BIOMASS": (6.0, 6.0)} for condition in objectives}
    result = compare_independent_and_joint_sparse_fba(
        parallel_pathway_model(),
        objectives=objectives,
        reaction_bounds=bounds,
        independent_lambda=0.1,
        joint_lambda=10.0,
    )

    independent_routes = {
        solution.condition: set(solution.active_by_flux) & {"ROUTE_A", "ROUTE_B"}
        for solution in result.independent
    }
    joint_routes = [
        set(solution.active_by_flux) & {"ROUTE_A", "ROUTE_B"}
        for solution in result.joint
    ]
    assert independent_routes == {
        "sample_prefers_A": {"ROUTE_A"},
        "sample_prefers_B": {"ROUTE_B"},
    }
    assert len({tuple(sorted(routes)) for routes in joint_routes}) == 1
    assert len(joint_routes[0]) == 1
    assert len(result.joint_active_union) < len(result.independent_active_union)


def test_joint_only_sparse_fba_preserves_condition_order() -> None:
    objectives = {
        "condition_A": {"ROUTE_B": 1.0},
        "condition_B": {"ROUTE_A": 1.0},
    }
    bounds = {condition: {"BIOMASS": (6.0, 6.0)} for condition in objectives}
    result = solve_joint_sparse_fba(
        parallel_pathway_model(),
        objectives=objectives,
        reaction_bounds=bounds,
        joint_lambda=10.0,
    )
    assert result.conditions == tuple(objectives)
    assert tuple(item.condition for item in result.joint) == tuple(objectives)
    assert result.to_dict()["joint_active_union_size"] == len(result.joint_active_union)
