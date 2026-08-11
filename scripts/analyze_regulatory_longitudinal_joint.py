#!/usr/bin/env python3
"""Summarize response-blind edge changes across all available joint receipts.

This is intentionally a descriptive, fail-closed post-processing step.  It
merges the four true-joint cohort receipts at one fixed lambda, aligns them by
run accession, and compares the primary runs listed in the Taylor family
design.  It does not use phenotype, treatment labels, or a solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


class LongitudinalError(ValueError):
    """Raised when a receipt or family design violates the input contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _edge_key(edge: Any, where: str) -> tuple[str, str, int]:
    if not isinstance(edge, dict):
        raise LongitudinalError(f"{where} must be an object")
    source, target, sign = edge.get("source"), edge.get("target"), edge.get("sign")
    if not isinstance(source, str) or not source or not isinstance(target, str) or not target:
        raise LongitudinalError(f"{where} has invalid source/target")
    if sign not in (-1, 1):
        raise LongitudinalError(f"{where} sign must be -1 or 1")
    return source, target, int(sign)


def _read_receipt(path: Path, expected_lambda: float) -> dict[str, dict[str, Any]]:
    try:
        root = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise LongitudinalError(f"cannot read receipt {path}: {error}") from error
    if not isinstance(root, dict) or root.get("status") != "completed":
        raise LongitudinalError(f"receipt {path} is not completed")
    method = root.get("method") if isinstance(root.get("method"), dict) else {}
    observed = method.get("lambda_nominal", method.get("lambda_reg_reported"))
    if observed is not None and abs(float(observed) - expected_lambda) > 1e-12:
        raise LongitudinalError(f"receipt {path} lambda {observed!r} != {expected_lambda}")
    rows = root.get("samples")
    if rows is None:
        rows = root.get("conditions")
    if not isinstance(rows, list) or not rows:
        raise LongitudinalError(f"receipt {path} has no samples/conditions")
    result: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise LongitudinalError(f"{path} record {index} is not an object")
        run = row.get("run_accession")
        if not isinstance(run, str) or not run or run in result:
            raise LongitudinalError(f"{path} record {index} has duplicate/invalid run")
        raw_edges = row.get("selected_edges")
        if not isinstance(raw_edges, list):
            raise LongitudinalError(f"{path} record {run} has no selected_edges")
        edges = {_edge_key(edge, f"{path}:{run}") for edge in raw_edges}
        result[run] = {
            "study_accession": row.get("study_accession", root.get("study_accession")),
            "status": str(row.get("status", "unknown")),
            "edges": edges,
        }
    return result


def _runs(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def analyze(receipts: list[Path], families: Path, output: Path, expected_lambda: float) -> dict[str, Any]:
    merged: dict[str, dict[str, Any]] = {}
    receipt_meta = []
    for path in receipts:
        records = _read_receipt(path, expected_lambda)
        overlap = set(merged).intersection(records)
        if overlap:
            raise LongitudinalError(f"duplicate runs across receipts: {sorted(overlap)[:3]}")
        merged.update(records)
        receipt_meta.append({"path": str(path), "sha256": _sha256(path), "record_count": len(records)})

    with families.open(encoding="utf-8", newline="") as handle:
        family_rows = list(csv.DictReader(handle, delimiter="\t"))
    pairs: list[dict[str, Any]] = []
    family_counts: Counter[str] = Counter()
    for family in family_rows:
        runs = _runs(family.get("primary_run_accessions", ""))
        available = [run for run in runs if run in merged]
        if len(available) < 2:
            continue
        baseline = available[0]
        baseline_edges = merged[baseline]["edges"]
        for current in available[1:]:
            current_edges = merged[current]["edges"]
            union = baseline_edges | current_edges
            pairs.append(
                {
                    "family_id": family.get("family_id", family.get("patient_id", "")),
                    "patient_id": family.get("patient_id", ""),
                    "analysis_role": family.get("analysis_role", ""),
                    "pair_status": family.get("pair_status", ""),
                    "baseline_run": baseline,
                    "current_run": current,
                    "baseline_study": merged[baseline]["study_accession"],
                    "current_study": merged[current]["study_accession"],
                    "baseline_status": merged[baseline]["status"],
                    "current_status": merged[current]["status"],
                    "baseline_edge_count": len(baseline_edges),
                    "current_edge_count": len(current_edges),
                    "gained_edge_count": len(current_edges - baseline_edges),
                    "lost_edge_count": len(baseline_edges - current_edges),
                    "edge_jaccard": len(current_edges & baseline_edges) / len(union) if union else 1.0,
                }
            )
            family_counts[str(family.get("family_id", ""))] += 1
    result = {
        "schema_version": "regulatory_longitudinal_joint.v1",
        "status": "completed",
        "response_blind": True,
        "lambda_nominal": expected_lambda,
        "receipt_count": len(receipts),
        "receipt_records": receipt_meta,
        "merged_run_count": len(merged),
        "family_design": {"path": str(families), "sha256": _sha256(families)},
        "pair_count": len(pairs),
        "family_pair_counts": dict(sorted(family_counts.items())),
        "status_counts": dict(sorted(Counter(f"{p['baseline_status']}__{p['current_status']}" for p in pairs).items())),
        "pairs": pairs,
        "claim_limit": "response-blind within-family descriptive rewiring only; no treatment causality or phenotype association",
        "completed_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise LongitudinalError(f"refusing to overwrite {output}")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", action="append", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lambda-nominal", type=float, default=0.001)
    args = parser.parse_args()
    result = analyze(args.receipt, args.families, args.output, args.lambda_nominal)
    print(json.dumps({"status": result["status"], "pair_count": result["pair_count"], "merged_run_count": result["merged_run_count"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
