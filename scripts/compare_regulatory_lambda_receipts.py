#!/usr/bin/env python3
"""Compare full-cohort response-blind regulatory receipts across lambda values."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any


def _jaccard(left: set[tuple[str, str, int]], right: set[tuple[str, str, int]]) -> float:
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def _edges(sample: dict[str, Any]) -> set[tuple[str, str, int]]:
    return {
        (str(edge.get("source")), str(edge.get("target")), int(edge.get("sign", 0)))
        for edge in sample.get("selected_edges", [])
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="append", required=True, metavar="LAMBDA=JSON")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipts: dict[str, dict[str, Any]] = {}
    for spec in args.pilot:
        if "=" not in spec:
            raise SystemExit("--pilot must be LAMBDA=JSON")
        label, raw_path = spec.split("=", 1)
        data = json.loads(Path(raw_path).read_text(encoding="utf-8"))
        if data.get("status") != "completed":
            raise SystemExit(f"{label}: receipt status is {data.get('status')!r}")
        samples = data.get("samples")
        if not isinstance(samples, list) or not samples:
            raise SystemExit(f"{label}: samples missing")
        receipts[label] = {
            "path": raw_path,
            "samples": samples,
            "sample_count": len(samples),
            "status_counts": dict(Counter(str(sample.get("status")) for sample in samples)),
            "edge_sets": [_edges(sample) for sample in samples],
        }
    labels = sorted(receipts, key=lambda value: float(value))
    baseline_label = "0" if "0" in receipts else labels[0]
    baseline = receipts[baseline_label]
    summary: dict[str, Any] = {}
    for label in labels:
        value = receipts[label]
        union = set().union(*value["edge_sets"])
        baseline_union = set().union(*baseline["edge_sets"])
        per_sample = [
            _jaccard(left, right)
            for left, right in zip(baseline["edge_sets"], value["edge_sets"])
        ] if len(value["edge_sets"]) == len(baseline["edge_sets"]) else []
        summary[label] = {
            "path": value["path"],
            "sample_count": value["sample_count"],
            "status_counts": value["status_counts"],
            "edge_union_size": len(union),
            "edge_jaccard_vs_baseline": _jaccard(baseline_union, union),
            "mean_sample_edge_jaccard_vs_baseline": (sum(per_sample) / len(per_sample)) if per_sample else None,
            "sample_edge_jaccards_vs_baseline": per_sample,
        }
    edge_frequency = Counter(
        edge for value in receipts.values() for edge in set().union(*value["edge_sets"])
    )
    result = {
        "schema_version": "regulatory_lambda_robustness.v1",
        "status": "completed",
        "response_blind": True,
        "baseline_lambda": baseline_label,
        "lambda_values": [float(label) for label in labels],
        "receipts": summary,
        "edge_frequency_across_lambda": [
            {"source": edge[0], "target": edge[1], "sign": edge[2], "lambda_count": count}
            for edge, count in edge_frequency.most_common()
        ],
        "claim_limit": "parameter-stability description only; no drug-response or causal interpretation",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "lambda_count": len(labels), "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
