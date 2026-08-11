from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "summarize_regulatory_lambda_all60.py"
SPEC = importlib.util.spec_from_file_location("regulatory_lambda_summary", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _receipt(lambda_value: float, shifted: bool = False):
    edge_a = ("A", "B", 1)
    edge_b = ("B", "C", -1)
    return {
        "path": f"lambda-{lambda_value}.json",
        "receipt_sha256": str(lambda_value),
        "source_sha256": {"expression": "same"},
        "samples": {
            "RUN1": {"status": "optimal", "edges": {edge_b if shifted else edge_a}},
            "RUN2": {"status": "blocked_no_selected_edges", "edges": set()},
        },
    }


def test_summary_labels_independent_solves_and_aligns_by_run():
    result = MODULE.summarize(
        {
            "STUDY-A": {0.0: _receipt(0.0), 0.1: _receipt(0.1, shifted=True)},
            "STUDY-B": {0.0: _receipt(0.0), 0.1: _receipt(0.1)},
        }
    )
    assert result["primary_sample_count"] == 4
    assert result["inference_scope"]["joint_multi_sample_inference"] is False
    assert result["cohorts"]["STUDY-A"]["lambda_summaries"]["0.1"][
        "mean_sample_edge_jaccard_vs_lambda0"
    ] == 0.5
    assert result["lambda_selection"]["status"] == "not_selected"


def test_summary_rejects_provenance_drift():
    changed = _receipt(0.1)
    changed["source_sha256"] = {"expression": "changed"}
    with pytest.raises(MODULE.ReceiptError, match="source_sha256 differs"):
        MODULE.summarize({"STUDY-A": {0.0: _receipt(0.0), 0.1: changed}})


def test_summary_rejects_sample_set_drift():
    changed = _receipt(0.1)
    changed["samples"].pop("RUN2")
    with pytest.raises(MODULE.ReceiptError, match="sample set differs"):
        MODULE.summarize({"STUDY-A": {0.0: _receipt(0.0), 0.1: changed}})
