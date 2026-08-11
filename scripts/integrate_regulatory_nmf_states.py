#!/usr/bin/env python3
"""Response-blind integration of regulatory CORNETO edges with NMF states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from scipy.stats import fisher_exact


class IntegrationError(ValueError):
    """Raised when an input fails the integration contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _edge_key(edge: Any, where: str) -> tuple[str, str, int]:
    if not isinstance(edge, dict):
        raise IntegrationError(f"{where} must be an object")
    source, target, sign = edge.get("source"), edge.get("target"), edge.get("sign")
    if not isinstance(source, str) or not source or not isinstance(target, str) or not target:
        raise IntegrationError(f"{where} has an invalid source or target")
    if sign not in (-1, 1):
        raise IntegrationError(f"{where} sign must be -1 or 1")
    return source, target, int(sign)


def _read_regulatory(path: Path, study: str) -> dict[str, dict[str, Any]]:
    try:
        root = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise IntegrationError(f"cannot read regulatory receipt: {error}") from error
    if not isinstance(root, dict) or root.get("status") != "completed":
        raise IntegrationError("regulatory receipt must have status='completed'")
    receipt_study = root.get("study", root.get("study_accession"))
    if receipt_study != study:
        raise IntegrationError(f"regulatory study {receipt_study!r} != {study!r}")
    # Legacy per-sample receipts call these records ``samples``.  The
    # normalized true-joint runner uses ``conditions`` for the same
    # run/edge/status records.
    samples = root.get("samples")
    if samples is None:
        samples = root.get("conditions")
    if not isinstance(samples, list) or not samples:
        raise IntegrationError("regulatory receipt has no samples/conditions")
    result: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(samples):
        if not isinstance(row, dict):
            raise IntegrationError(f"samples[{index}] must be an object")
        run = row.get("run_accession")
        status = row.get("status")
        if not isinstance(run, str) or not run or run in result:
            raise IntegrationError(f"samples[{index}] has invalid/duplicate run_accession")
        if not isinstance(status, str) or not status:
            raise IntegrationError(f"samples[{index}] has no status")
        edges_raw = row.get("selected_edges")
        if not isinstance(edges_raw, list):
            raise IntegrationError(f"samples[{index}].selected_edges must be an array")
        edges = {_edge_key(edge, f"samples[{index}].selected_edges") for edge in edges_raw}
        if len(edges) != len(edges_raw):
            raise IntegrationError(f"samples[{index}] has duplicate selected edges")
        result[run] = {"status": status, "edges": edges}
    return result


def _read_assignments(path: Path, study: str) -> dict[str, dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except OSError as error:
        raise IntegrationError(f"cannot read NMF assignments: {error}") from error
    state_column = "independent_state" if rows and "independent_state" in rows[0] else "state"
    required = {"study_accession", "run_accession", state_column}
    if not rows or not required.issubset(rows[0]):
        raise IntegrationError("NMF assignments are empty or missing required columns")
    result: dict[str, dict[str, str]] = {}
    for index, row in enumerate(rows, start=2):
        run, state = row["run_accession"], row[state_column]
        if row["study_accession"] != study:
            raise IntegrationError(f"assignments line {index} has wrong study")
        if not run or not state or run in result:
            raise IntegrationError(f"assignments line {index} has invalid/duplicate identifiers")
        result[run] = {"state": state, **row}
    return result


def _jaccard(left: set[tuple[str, str, int]], right: set[tuple[str, str, int]]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def _bh(rows: list[dict[str, Any]]) -> None:
    order = sorted(range(len(rows)), key=lambda index: (rows[index]["p_value"], index))
    adjusted = [1.0] * len(rows)
    running = 1.0
    for rank_index in range(len(order) - 1, -1, -1):
        index = order[rank_index]
        running = min(running, rows[index]["p_value"] * len(rows) / (rank_index + 1))
        adjusted[index] = running
    for row, value in zip(rows, adjusted, strict=True):
        row["q_value_bh"] = value


def integrate(regulatory_path: Path, assignments_path: Path, study: str) -> dict[str, Any]:
    regulatory = _read_regulatory(regulatory_path, study)
    assignments = _read_assignments(assignments_path, study)
    matched = sorted(set(regulatory) & set(assignments))
    if len(matched) < 2:
        raise IntegrationError("fewer than two samples overlap")
    states = sorted({assignments[run]["state"] for run in matched})
    if len(states) < 2:
        raise IntegrationError("fewer than two NMF states overlap")
    runs_by_state = {
        state: [run for run in matched if assignments[run]["state"] == state]
        for state in states
    }
    unions = {
        state: set().union(*(regulatory[run]["edges"] for run in runs))
        for state, runs in runs_by_state.items()
    }
    summaries = []
    for state, runs in runs_by_state.items():
        burdens = [len(regulatory[run]["edges"]) for run in runs]
        summaries.append(
            {
                "state": state,
                "sample_count": len(runs),
                "status_counts": dict(
                    sorted(Counter(regulatory[run]["status"] for run in runs).items())
                ),
                "edge_burden_median": statistics.median(burdens),
                "edge_burden_min": min(burdens),
                "edge_burden_max": max(burdens),
                "edge_union_count": len(unions[state]),
            }
        )
    pairwise = []
    for index, left in enumerate(states):
        for right in states[index + 1 :]:
            pairwise.append(
                {
                    "left": left,
                    "right": right,
                    "edge_intersection": len(unions[left] & unions[right]),
                    "edge_union": len(unions[left] | unions[right]),
                    "edge_jaccard": _jaccard(unions[left], unions[right]),
                }
            )
    all_edges = sorted(set().union(*(regulatory[run]["edges"] for run in matched)))
    tests: list[dict[str, Any]] = []
    for edge in all_edges:
        for state in states:
            inside = runs_by_state[state]
            inside_set = set(inside)
            outside = [run for run in matched if run not in inside_set]
            a = sum(edge in regulatory[run]["edges"] for run in inside)
            c = sum(edge in regulatory[run]["edges"] for run in outside)
            _, p_value = fisher_exact([[a, len(inside) - a], [c, len(outside) - c]])
            tests.append(
                {
                    "source": edge[0],
                    "target": edge[1],
                    "sign": edge[2],
                    "state": state,
                    "state_present": a,
                    "state_total": len(inside),
                    "other_present": c,
                    "other_total": len(outside),
                    "prevalence_difference": a / len(inside) - c / len(outside),
                    "p_value": float(p_value),
                }
            )
    _bh(tests)
    tests.sort(
        key=lambda row: (
            row["q_value_bh"],
            row["p_value"],
            -abs(row["prevalence_difference"]),
        )
    )
    return {
        "schema_version": "regulatory_nmf_integration.v1",
        "status": "completed",
        "response_blind": True,
        "study_accession": study,
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "inputs": {
            "regulatory_receipt": str(regulatory_path),
            "regulatory_sha256": _sha256(regulatory_path),
            "nmf_assignments": str(assignments_path),
            "nmf_assignments_sha256": _sha256(assignments_path),
        },
        "coverage": {
            "regulatory_samples": len(regulatory),
            "nmf_samples": len(assignments),
            "matched_samples": len(matched),
            "regulatory_only": sorted(set(regulatory) - set(assignments)),
            "nmf_only": sorted(set(assignments) - set(regulatory)),
        },
        "state_summaries": summaries,
        "pairwise_state_edge_overlap": pairwise,
        "edge_state_tests": tests,
        "multiple_testing": (
            "two-sided Fisher exact tests; Benjamini-Hochberg across all edge-state tests"
        ),
        "claim_limit": (
            "response-blind descriptive association between independently derived NMF states and "
            "CORNETO-selected regulatory edges; no Taxol-response, causal, or official "
            "Barnes-subtype claim"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regulatory", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = integrate(args.regulatory, args.assignments, args.study)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    except (IntegrationError, OSError) as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": "completed",
                "matched_samples": result["coverage"]["matched_samples"],
                "states": len(result["state_summaries"]),
                "edge_state_tests": len(result["edge_state_tests"]),
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
