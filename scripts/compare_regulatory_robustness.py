#!/usr/bin/env python3
"""Fail-closed comparison of response-blind regulatory sensitivity receipts."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


class RobustnessError(ValueError):
    """Raised when receipts cannot be compared fairly."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise RobustnessError(f"receipt must be LABEL=PATH, got {spec!r}")
    label, raw_path = spec.split("=", 1)
    if not label or not raw_path:
        raise RobustnessError(f"receipt must be LABEL=PATH, got {spec!r}")
    return label, Path(raw_path)


def _edge(edge: Any, where: str) -> tuple[str, str, int]:
    if not isinstance(edge, dict):
        raise RobustnessError(f"{where} must be an object")
    source, target, sign = edge.get("source"), edge.get("target"), edge.get("sign")
    if not isinstance(source, str) or not source or not isinstance(target, str) or not target:
        raise RobustnessError(f"{where} has invalid nodes")
    if sign not in (-1, 1):
        raise RobustnessError(f"{where} has invalid sign")
    return source, target, int(sign)


def _read(label: str, path: Path, study: str) -> dict[str, Any]:
    try:
        root = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RobustnessError(f"{label}: cannot read receipt: {error}") from error
    if not isinstance(root, dict) or root.get("status") != "completed":
        raise RobustnessError(f"{label}: receipt status must be 'completed'")
    receipt_study = root.get("study", root.get("study_accession"))
    if receipt_study != study or root.get("method", {}).get("response_blind") is not True:
        raise RobustnessError(f"{label}: study/response-blind contract failed")
    samples_raw = root.get("samples")
    if not isinstance(samples_raw, list) or not samples_raw:
        raise RobustnessError(f"{label}: samples must be a non-empty array")
    samples: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(samples_raw):
        if not isinstance(row, dict):
            raise RobustnessError(f"{label}: sample {index} is not an object")
        run, status = row.get("run_accession"), row.get("status")
        if not isinstance(run, str) or not run or run in samples:
            raise RobustnessError(f"{label}: invalid/duplicate run_accession")
        if not isinstance(status, str) or not status:
            raise RobustnessError(f"{label}: sample {run} has no status")
        raw_edges = row.get("selected_edges")
        if not isinstance(raw_edges, list):
            raise RobustnessError(f"{label}: sample {run} selected_edges is not an array")
        edges = {_edge(item, f"{label}.{run}.selected_edges") for item in raw_edges}
        if len(edges) != len(raw_edges):
            raise RobustnessError(f"{label}: sample {run} has duplicate edges")
        samples[run] = {"status": status, "edges": edges}
    return {
        "label": label,
        "path": str(path),
        "sha256": _sha256(path),
        "method": root.get("method"),
        "source_sha256": root.get("source_sha256"),
        "samples": samples,
    }


def _jaccard(left: set[tuple[str, str, int]], right: set[tuple[str, str, int]]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def compare(specs: list[str], study: str) -> dict[str, Any]:
    parsed = [_parse(spec) for spec in specs]
    if len(parsed) < 2 or len({label for label, _ in parsed}) != len(parsed):
        raise RobustnessError("at least two uniquely labelled receipts are required")
    receipts = [_read(label, path, study) for label, path in parsed]
    reference_runs = set(receipts[0]["samples"])
    reference_sources = receipts[0]["source_sha256"]
    for receipt in receipts[1:]:
        if set(receipt["samples"]) != reference_runs:
            raise RobustnessError(f"{receipt['label']}: sample set differs from reference")
        if receipt["source_sha256"] != reference_sources:
            raise RobustnessError(f"{receipt['label']}: source checksums differ from reference")

    summaries = []
    union_by_label: dict[str, set[tuple[str, str, int]]] = {}
    edge_policy_frequency: Counter[tuple[str, str, int]] = Counter()
    for receipt in receipts:
        samples = receipt["samples"]
        union = set().union(*(row["edges"] for row in samples.values()))
        union_by_label[receipt["label"]] = union
        edge_policy_frequency.update(union)
        statuses = Counter(row["status"] for row in samples.values())
        edge_counts = Counter(edge for row in samples.values() for edge in row["edges"])
        summaries.append(
            {
                "label": receipt["label"],
                "path": receipt["path"],
                "sha256": receipt["sha256"],
                "method": receipt["method"],
                "sample_count": len(samples),
                "status_counts": dict(sorted(statuses.items())),
                "blocked_count": sum(count for status, count in statuses.items() if status.startswith("blocked")),
                "edge_union_count": len(union),
                "edge_selection_frequency": [
                    {
                        "source": edge[0],
                        "target": edge[1],
                        "sign": edge[2],
                        "sample_count": count,
                    }
                    for edge, count in sorted(edge_counts.items(), key=lambda item: (-item[1], item[0]))
                ],
            }
        )

    pairwise = []
    for index, left in enumerate(receipts):
        for right in receipts[index + 1 :]:
            per_sample = []
            for run in sorted(reference_runs):
                left_row, right_row = left["samples"][run], right["samples"][run]
                per_sample.append(
                    {
                        "run_accession": run,
                        "left_status": left_row["status"],
                        "right_status": right_row["status"],
                        "left_edge_count": len(left_row["edges"]),
                        "right_edge_count": len(right_row["edges"]),
                        "edge_jaccard": _jaccard(left_row["edges"], right_row["edges"]),
                    }
                )
            jaccards = [row["edge_jaccard"] for row in per_sample]
            left_union = union_by_label[left["label"]]
            right_union = union_by_label[right["label"]]
            pairwise.append(
                {
                    "left": left["label"],
                    "right": right["label"],
                    "status_agreement_count": sum(
                        row["left_status"] == row["right_status"] for row in per_sample
                    ),
                    "sample_count": len(per_sample),
                    "per_sample_edge_jaccard_median": statistics.median(jaccards),
                    "per_sample_edge_jaccard_min": min(jaccards),
                    "per_sample_edge_jaccard_max": max(jaccards),
                    "pooled_edge_jaccard": _jaccard(left_union, right_union),
                    "per_sample": per_sample,
                }
            )
    return {
        "schema_version": "regulatory_robustness.v1",
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "study_accession": study,
        "response_blind": True,
        "sample_count": len(reference_runs),
        "policies": summaries,
        "pairwise": pairwise,
        "edge_policy_frequency": [
            {
                "source": edge[0],
                "target": edge[1],
                "sign": edge[2],
                "policy_count": count,
            }
            for edge, count in sorted(
                edge_policy_frequency.items(), key=lambda item: (-item[1], item[0])
            )
        ],
        "claim_limit": (
            "response-blind technical sensitivity to regularisation policy; no phenotype, "
            "treatment-response, or causal claim"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", action="append", required=True, metavar="LABEL=PATH")
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compare(args.receipt, args.study)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    except (OSError, RobustnessError) as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": "completed",
                "sample_count": result["sample_count"],
                "policy_count": len(result["policies"]),
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
