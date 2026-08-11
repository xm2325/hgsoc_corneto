import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "summarize_regulatory_alternative_optima.py"
SPEC = importlib.util.spec_from_file_location("summarize_regulatory_alternative_optima", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def _receipt(run: str, edges: int) -> dict:
    return {
        "schema_version": "regulatory_alternative_optima.v1",
        "study": "E-MTAB-X",
        "method": {"lambda_reg": 0.1},
        "source_sha256": {"expression": "abc"},
        "samples": [
            {
                "run_accession": run,
                "status": "completed",
                "solution_count": 3,
                "accepted_alternative_count": 2,
                "incumbent_edge_count": edges,
                "edge_union_count": edges,
                "edge_intersection_count": edges,
                "core_edge_count": edges,
                "mean_pairwise_jaccard": 1.0,
                "min_pairwise_jaccard": 1.0,
                "mean_incumbent_jaccard": 1.0,
                "mean_edge_entropy_bits": 0.0,
            }
        ],
    }


def test_summary_reclassifies_zero_edge_incumbent(tmp_path):
    paths = []
    for run, edges in (("RUN1", 2), ("RUN2", 0)):
        path = tmp_path / f"{run}.json"
        path.write_text(json.dumps(_receipt(run, edges)), encoding="utf-8")
        paths.append(path)

    result = MODULE.summarize(paths, expected_samples=2, study="E-MTAB-X")

    assert result["counts"]["usable_nonempty_samples"] == 1
    assert result["counts"]["status_counts"] == {
        "blocked_no_selected_edges": 1,
        "completed": 1,
    }
    assert result["ensemble_metrics_nonempty_samples"]["median_incumbent_edge_count"] == 2


def test_summary_prefers_successful_retry(tmp_path):
    failed = _receipt("RUN1", 0)
    failed["samples"][0]["status"] = "error"
    original = tmp_path / "original.json"
    original.write_text(json.dumps(failed), encoding="utf-8")
    retry = tmp_path / "retry.json"
    retry.write_text(json.dumps(_receipt("RUN1", 4)), encoding="utf-8")

    result = MODULE.summarize([original, retry], expected_samples=1, study="E-MTAB-X")

    assert result["counts"]["input_receipts"] == 2
    assert result["counts"]["superseded_receipts"] == 1
    assert result["samples"][0]["status"] == "completed"
    assert result["samples"][0]["incumbent_edge_count"] == 4
