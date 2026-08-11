import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_regulatory_alternative_optima.py"
SPEC = importlib.util.spec_from_file_location("run_regulatory_alternative_optima", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_summarize_edge_ensemble_counts_and_jaccard():
    result = MODULE.summarize_edge_ensemble(
        [{"A", "B"}, {"A", "C"}, {"A", "B"}], core_frequency=2 / 3
    )

    assert result["solution_count"] == 3
    assert result["accepted_alternative_count"] == 2
    assert result["incumbent_edge_count"] == 2
    assert result["edge_union_count"] == 3
    assert result["edge_intersection_count"] == 1
    assert result["core_edge_count"] == 2
    assert result["mean_pairwise_jaccard"] == (1 / 3 + 1 + 1 / 3) / 3
    assert [item["edge"] for item in result["edge_frequencies"]] == ["A", "B", "C"]


def test_empty_ensemble_is_rejected():
    try:
        MODULE.summarize_edge_ensemble([], core_frequency=0.8)
    except ValueError as error:
        assert "empty" in str(error)
    else:
        raise AssertionError("empty ensemble should fail")
