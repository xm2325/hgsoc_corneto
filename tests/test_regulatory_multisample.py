from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

SCRIPT = Path(__file__).parents[1] / "scripts" / "run_corneto_regulatory_multisample.py"
SPEC = importlib.util.spec_from_file_location("run_corneto_regulatory_multisample", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_effective_lambda_normalizes_condition_fit_sum() -> None:
    assert MODULE._effective_lambda(0.25, 60, "mean_fit") == 15.0
    assert MODULE._effective_lambda(0.25, 9, "mean_fit") == 2.25
    assert MODULE._effective_lambda(0.25, 60, "raw") == 0.25
    with pytest.raises(ValueError, match="non-negative"):
        MODULE._effective_lambda(-0.1, 60, "mean_fit")


def test_edge_records_preserve_condition_columns() -> None:
    class Graph:
        E = [({"A"}, {"B"}), ({"B"}, {"C"})]

        @staticmethod
        def get_attr_edge(index: int) -> dict[str, int]:
            return {"interaction": 1 if index == 0 else -1}

    method = SimpleNamespace(processed_graph=Graph())
    selected = MODULE._edge_records(
        method,
        np.asarray([[1, 0], [1, 1]], dtype=float),
        ["condition_a", "condition_b"],
    )
    assert selected["condition_a"] == [
        {"source": "A", "target": "B", "sign": 1},
        {"source": "B", "target": "C", "sign": -1},
    ]
    assert selected["condition_b"] == [
        {"source": "B", "target": "C", "sign": -1}
    ]
