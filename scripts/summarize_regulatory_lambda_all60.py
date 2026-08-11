#!/usr/bin/env python3
"""Audit and summarize response-blind regulatory lambda receipts.

The current regulatory runner solves one CARNIVAL problem per sample.  This
script therefore calls the 60-OCM result a *pooled descriptive aggregation*,
not joint multi-sample inference.  It validates sample identity and provenance
before comparing edge stability across lambda values and cohorts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


Edge = tuple[str, str, int]


class ReceiptError(ValueError):
    """Raised when a regulatory receipt violates the comparison contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _jaccard(left: set[Edge], right: set[Edge]) -> float:
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def _mean(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def _median(values: list[float | int]) -> float | None:
    return statistics.median(values) if values else None


def _edge(raw: Any, where: str) -> Edge:
    if not isinstance(raw, dict):
        raise ReceiptError(f"{where} must be an object")
    source, target, sign = raw.get("source"), raw.get("target"), raw.get("sign")
    if not isinstance(source, str) or not source or not isinstance(target, str) or not target:
        raise ReceiptError(f"{where} has an invalid source or target")
    if sign not in (-1, 1):
        raise ReceiptError(f"{where} sign must be -1 or 1")
    return source, target, int(sign)


def _load_receipt(path: Path, expected_study: str, expected_lambda: float) -> dict[str, Any]:
    try:
        root = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReceiptError(f"cannot read {path}: {error}") from error
    if not isinstance(root, dict) or root.get("status") not in {"completed", "blocked"}:
        raise ReceiptError(f"{path}: receipt status must be 'completed' or 'blocked'")
    if root.get("study") != expected_study:
        raise ReceiptError(f"{path}: study {root.get('study')!r} != {expected_study!r}")
    if root.get("primary_only") is not True:
        raise ReceiptError(f"{path}: primary_only must be true")
    method = root.get("method")
    if not isinstance(method, dict) or method.get("response_blind") is not True:
        raise ReceiptError(f"{path}: response_blind method flag missing")
    recorded_lambda = method.get("lambda_reg")
    if not isinstance(recorded_lambda, (int, float)) or not math.isclose(
        float(recorded_lambda), expected_lambda, rel_tol=0.0, abs_tol=1e-12
    ):
        raise ReceiptError(
            f"{path}: method.lambda_reg={recorded_lambda!r} != {expected_lambda!r}"
        )
    raw_samples = root.get("samples")
    if not isinstance(raw_samples, list) or not raw_samples:
        raise ReceiptError(f"{path}: samples missing")
    samples: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(raw_samples):
        if not isinstance(row, dict):
            raise ReceiptError(f"{path}: samples[{index}] must be an object")
        run = row.get("run_accession")
        if not isinstance(run, str) or not run or run in samples:
            raise ReceiptError(f"{path}: samples[{index}] has invalid/duplicate run_accession")
        status = row.get("status")
        if not isinstance(status, str) or not status:
            raise ReceiptError(f"{path}: samples[{index}] has no status")
        raw_edges = row.get("selected_edges")
        if not isinstance(raw_edges, list):
            raise ReceiptError(f"{path}: samples[{index}].selected_edges must be an array")
        edges = {_edge(value, f"{path}: samples[{index}].selected_edges") for value in raw_edges}
        if len(edges) != len(raw_edges):
            raise ReceiptError(f"{path}: samples[{index}] has duplicate selected edges")
        samples[run] = {"status": status, "edges": edges}
    sample_statuses = {row["status"] for row in samples.values()}
    if "error" in sample_statuses:
        raise ReceiptError(f"{path}: one or more sample solves have status='error'")
    if root["status"] == "blocked" and any(
        status in {"optimal", "optimal_inaccurate"} for status in sample_statuses
    ):
        raise ReceiptError(f"{path}: blocked receipt contains an optimal sample")
    source_sha = root.get("source_sha256")
    if not isinstance(source_sha, dict) or not source_sha:
        raise ReceiptError(f"{path}: source_sha256 missing")
    return {
        "path": str(path),
        "receipt_sha256": _sha256(path),
        "source_sha256": source_sha,
        "receipt_status": root["status"],
        "samples": samples,
    }


def _study_summary(values: dict[float, dict[str, Any]]) -> dict[str, Any]:
    if 0.0 not in values:
        raise ReceiptError("lambda=0 baseline receipt missing")
    lambdas = sorted(values)
    baseline = values[0.0]
    runs = set(baseline["samples"])
    source_sha = baseline["source_sha256"]
    for value in values.values():
        if set(value["samples"]) != runs:
            missing = sorted(runs - set(value["samples"]))
            extra = sorted(set(value["samples"]) - runs)
            raise ReceiptError(f"sample set differs from lambda=0; missing={missing}, extra={extra}")
        if value["source_sha256"] != source_sha:
            raise ReceiptError("source_sha256 differs across lambda receipts")

    rows: dict[str, Any] = {}
    previous_union: set[Edge] | None = None
    for lambda_value in lambdas:
        value = values[lambda_value]
        samples = value["samples"]
        union = set().union(*(samples[run]["edges"] for run in sorted(runs)))
        baseline_union = set().union(
            *(baseline["samples"][run]["edges"] for run in sorted(runs))
        )
        per_sample = [
            _jaccard(baseline["samples"][run]["edges"], samples[run]["edges"])
            for run in sorted(runs)
        ]
        burdens = [len(samples[run]["edges"]) for run in sorted(runs)]
        status_transitions = Counter(
            f"{baseline['samples'][run]['status']} -> {samples[run]['status']}"
            for run in sorted(runs)
        )
        rows[format(lambda_value, ".12g")] = {
            "path": value["path"],
            "receipt_sha256": value["receipt_sha256"],
            "receipt_status": value["receipt_status"],
            "sample_count": len(runs),
            "status_counts": dict(sorted(Counter(samples[run]["status"] for run in runs).items())),
            "status_transitions_vs_lambda0": dict(sorted(status_transitions.items())),
            "edge_union_size": len(union),
            "edge_jaccard_vs_lambda0": _jaccard(baseline_union, union),
            "edge_jaccard_vs_previous_lambda": (
                _jaccard(previous_union, union) if previous_union is not None else None
            ),
            "sample_edge_burden_median": _median(burdens),
            "sample_edge_burden_min": min(burdens),
            "sample_edge_burden_max": max(burdens),
            "mean_sample_edge_jaccard_vs_lambda0": _mean(per_sample),
            "median_sample_edge_jaccard_vs_lambda0": _median(per_sample),
            "min_sample_edge_jaccard_vs_lambda0": min(per_sample),
        }
        previous_union = union

    sample_stability = []
    for run in sorted(runs):
        edge_sets = [values[value]["samples"][run]["edges"] for value in lambdas]
        pairwise = [
            _jaccard(edge_sets[left], edge_sets[right])
            for left in range(len(edge_sets))
            for right in range(left + 1, len(edge_sets))
        ]
        sample_stability.append(
            {
                "run_accession": run,
                "status_count_across_lambda": len(
                    {values[value]["samples"][run]["status"] for value in lambdas}
                ),
                "edge_intersection_size": len(set.intersection(*edge_sets)) if edge_sets else 0,
                "edge_union_size": len(set.union(*edge_sets)) if edge_sets else 0,
                "mean_pairwise_edge_jaccard": _mean(pairwise),
                "edge_burden_range": [min(map(len, edge_sets)), max(map(len, edge_sets))],
            }
        )

    union_frequency: Counter[Edge] = Counter()
    occurrence_frequency: Counter[Edge] = Counter()
    for lambda_value in lambdas:
        edge_sets = [values[lambda_value]["samples"][run]["edges"] for run in sorted(runs)]
        for edge in set().union(*edge_sets):
            union_frequency[edge] += 1
        for edge_set in edge_sets:
            occurrence_frequency.update(edge_set)
    edge_stability = [
        {
            "source": edge[0],
            "target": edge[1],
            "sign": edge[2],
            "lambda_union_count": union_frequency[edge],
            "lambda_union_fraction": union_frequency[edge] / len(lambdas),
            "sample_lambda_occurrence_count": occurrence_frequency[edge],
            "sample_lambda_occurrence_fraction": occurrence_frequency[edge]
            / (len(runs) * len(lambdas)),
        }
        for edge in sorted(
            union_frequency,
            key=lambda edge: (-union_frequency[edge], -occurrence_frequency[edge], edge),
        )
    ]
    return {
        "sample_count": len(runs),
        "lambda_values": lambdas,
        "source_sha256": source_sha,
        "lambda_summaries": rows,
        "sample_stability": sample_stability,
        "edge_stability": edge_stability,
    }


def summarize(receipts: dict[str, dict[float, dict[str, Any]]]) -> dict[str, Any]:
    studies = sorted(receipts)
    if not studies:
        raise ReceiptError("no receipts provided")
    expected_lambdas = sorted(receipts[studies[0]])
    for study in studies[1:]:
        if sorted(receipts[study]) != expected_lambdas:
            raise ReceiptError(f"{study}: lambda grid differs across cohorts")
    by_study = {study: _study_summary(receipts[study]) for study in studies}

    pooled_rows: dict[str, Any] = {}
    pairwise_rows: list[dict[str, Any]] = []
    baseline_samples = {
        (study, run): row
        for study in studies
        for run, row in receipts[study][0.0]["samples"].items()
    }
    for lambda_value in expected_lambdas:
        current_samples = {
            (study, run): row
            for study in studies
            for run, row in receipts[study][lambda_value]["samples"].items()
        }
        union = set().union(*(row["edges"] for row in current_samples.values()))
        baseline_union = set().union(*(row["edges"] for row in baseline_samples.values()))
        per_sample = [
            _jaccard(baseline_samples[key]["edges"], current_samples[key]["edges"])
            for key in sorted(current_samples)
        ]
        burdens = [len(row["edges"]) for row in current_samples.values()]
        blocked = sum(row["status"] == "blocked_no_selected_edges" for row in current_samples.values())
        pooled_rows[format(lambda_value, ".12g")] = {
            "sample_count": len(current_samples),
            "status_counts": dict(
                sorted(Counter(row["status"] for row in current_samples.values()).items())
            ),
            "blocked_no_selected_edges_fraction": blocked / len(current_samples),
            "edge_union_size": len(union),
            "edge_jaccard_vs_lambda0": _jaccard(baseline_union, union),
            "sample_edge_burden_median": _median(burdens),
            "mean_sample_edge_jaccard_vs_lambda0": _mean(per_sample),
            "median_sample_edge_jaccard_vs_lambda0": _median(per_sample),
        }
        cohort_unions = {
            study: set().union(
                *(row["edges"] for row in receipts[study][lambda_value]["samples"].values())
            )
            for study in studies
        }
        for left_index, left in enumerate(studies):
            for right in studies[left_index + 1 :]:
                pairwise_rows.append(
                    {
                        "lambda_reg": lambda_value,
                        "left": left,
                        "right": right,
                        "edge_jaccard": _jaccard(cohort_unions[left], cohort_unions[right]),
                        "edge_intersection_size": len(cohort_unions[left] & cohort_unions[right]),
                        "edge_union_size": len(cohort_unions[left] | cohort_unions[right]),
                    }
                )

    return {
        "schema_version": "regulatory_lambda_all_primary_audit.v2",
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "response_blind": True,
        "inference_scope": {
            "runner_behavior": "one independent CARNIVAL solve per OCM and lambda",
            "pooled_result": "descriptive aggregation of independent per-sample solves",
            "joint_multi_sample_inference": False,
            "important_limitation": (
                "these receipts do not test CORNETO joint multi-sample regularization; "
                "cohort preprocessing also changes expression z-scores and selected priors"
            ),
        },
        "study_count": len(studies),
        "primary_sample_count": sum(value["sample_count"] for value in by_study.values()),
        "lambda_values": expected_lambdas,
        "cohorts": by_study,
        "pooled_descriptive_aggregation": pooled_rows,
        "cohort_stratified_edge_overlap": pairwise_rows,
        "lambda_selection": {
            "status": "not_selected",
            "reason": (
                "edge stability alone is insufficient; predeclared fit/agreement criteria and "
                "alternative-solution sampling are required"
            ),
        },
        "claim_limit": (
            "response-blind parameter and edge-stability audit only; no joint multi-sample, "
            "drug-response, causal, or Barnes-subtype interpretation"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", action="append", required=True, metavar="STUDY|LAMBDA|JSON")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipts: dict[str, dict[float, dict[str, Any]]] = {}
    try:
        for spec in args.receipt:
            pieces = spec.split("|", 2)
            if len(pieces) != 3:
                raise ReceiptError("--receipt must be STUDY|LAMBDA|JSON")
            study, raw_lambda, raw_path = pieces
            lambda_value = float(raw_lambda)
            if lambda_value in receipts.setdefault(study, {}):
                raise ReceiptError(f"duplicate receipt for {study}, lambda={lambda_value}")
            receipts[study][lambda_value] = _load_receipt(
                Path(raw_path), study, lambda_value
            )
        result = summarize(receipts)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    except (ReceiptError, OSError, ValueError) as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": "completed",
                "study_count": result["study_count"],
                "primary_sample_count": result["primary_sample_count"],
                "lambda_count": len(result["lambda_values"]),
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
