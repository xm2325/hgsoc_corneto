import argparse
import csv
import json
from pathlib import Path

import pytest

from scripts.compare_pooled_stratified import ComparisonError, compare


def _tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _regulatory(path: Path, samples: dict[str, list[tuple[str, str, int]]], joint: bool) -> None:
    path.write_text(
        json.dumps(
            {
                "status": "completed",
                "method": {"lambda_reg": 0.1, "joint_multi_sample": joint},
                "samples": [
                    {
                        "run_accession": run,
                        "status": "optimal",
                        "selected_edges": [
                            {"source": source, "target": target, "sign": sign}
                            for source, target, sign in edges
                        ],
                    }
                    for run, edges in samples.items()
                ],
            }
        ),
        encoding="utf-8",
    )


def _args(tmp_path: Path, require_joint: bool = True) -> argparse.Namespace:
    manifest = tmp_path / "manifest.tsv"
    _tsv(
        manifest,
        [
            {"run_accession": "r1", "study_accession": "A", "canonical_ocm_id": "o1", "patient_id": "p1", "primary_cohort_eligible": "true"},
            {"run_accession": "r2", "study_accession": "A", "canonical_ocm_id": "o2", "patient_id": "p2", "primary_cohort_eligible": "true"},
            {"run_accession": "r3", "study_accession": "B", "canonical_ocm_id": "o3", "patient_id": "p3", "primary_cohort_eligible": "true"},
            {"run_accession": "r4", "study_accession": "B", "canonical_ocm_id": "o4", "patient_id": "p4", "primary_cohort_eligible": "true"},
        ],
    )
    pooled_reg = tmp_path / "pooled.json"
    cohort_a = tmp_path / "a.json"
    cohort_b = tmp_path / "b.json"
    edge_x = ("X", "Y", 1)
    edge_z = ("Z", "Y", -1)
    _regulatory(pooled_reg, {"r1": [edge_x], "r2": [edge_x], "r3": [edge_z], "r4": [edge_z]}, True)
    _regulatory(cohort_a, {"r1": [edge_x], "r2": [edge_z]}, False)
    _regulatory(cohort_b, {"r3": [edge_z], "r4": [edge_z]}, False)
    pooled_nmf = tmp_path / "pooled.tsv"
    nmf_a = tmp_path / "nmf_a.tsv"
    nmf_b = tmp_path / "nmf_b.tsv"
    _tsv(pooled_nmf, [{"run_accession": "r1", "pooled_state": "P1"}, {"run_accession": "r2", "pooled_state": "P2"}, {"run_accession": "r3", "pooled_state": "P1"}, {"run_accession": "r4", "pooled_state": "P2"}])
    _tsv(nmf_a, [{"run_accession": "r1", "independent_state": "A2"}, {"run_accession": "r2", "independent_state": "A1"}])
    _tsv(nmf_b, [{"run_accession": "r3", "independent_state": "B7"}, {"run_accession": "r4", "independent_state": "B8"}])
    return argparse.Namespace(
        manifest=manifest,
        pooled_regulatory=[f"0.1={pooled_reg}"],
        cohort_regulatory=[f"A|0.1|{cohort_a}", f"B|0.1|{cohort_b}"],
        require_joint_pooled=require_joint,
        pooled_nmf=pooled_nmf,
        cohort_nmf=[f"A|{nmf_a}", f"B|{nmf_b}"],
        pooled_state_column=None,
        cohort_state_column=None,
    )


def test_compare_edges_and_align_nmf_states(tmp_path: Path) -> None:
    result = compare(_args(tmp_path))
    row = result["regulatory"]["lambda_comparisons"][0]
    assert row["matched_sample_count"] == 4
    assert row["pooled_edge_union_size"] == 2
    assert row["cohort_merged_edge_union_size"] == 2
    assert row["edge_union_jaccard"] == 1.0
    assert row["mean_matched_sample_edge_jaccard"] == pytest.approx(0.75)
    assert result["manifest"]["primary_run_count"] == 4
    alignment = result["nmf"]["state_alignment"]
    assert len(alignment) == 4
    assert all(row["state_jaccard"] == 1.0 for row in alignment)
    assert all(row["adjusted_rand_index"] == 1.0 for row in result["nmf"]["study_summaries"])


def test_joint_contract_rejects_independent_pooled_receipt(tmp_path: Path) -> None:
    args = _args(tmp_path)
    pooled_path = Path(args.pooled_regulatory[0].split("=", 1)[1])
    data = json.loads(pooled_path.read_text(encoding="utf-8"))
    data["method"]["joint_multi_sample"] = False
    pooled_path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ComparisonError, match="joint_multi_sample"):
        compare(args)
