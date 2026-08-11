#!/usr/bin/env python3
"""Compare pooled-60 and patient-balanced true-joint regulatory receipts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path) -> dict[str, set[tuple[str, str, int]]]:
    root = json.loads(path.read_text(encoding="utf-8"))
    if root.get("status") != "completed":
        raise ValueError(f"{path} is not completed")
    rows = root.get("conditions")
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{path} has no conditions")
    return {
        row["run_accession"]: {(e["source"], e["target"], int(e["sign"])) for e in row.get("selected_edges", [])}
        for row in rows
    }


def _jaccard(left: set, right: set) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pooled", type=Path, required=True)
    parser.add_argument("--patient-balanced", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError(f"refusing to overwrite {args.output}")
    pooled, balanced = _read(args.pooled), _read(args.patient_balanced)
    common = sorted(set(pooled) & set(balanced))
    if len(common) != 52:
        raise ValueError(f"expected 52 common runs, found {len(common)}")
    pooled_union = set().union(*(pooled[run] for run in common))
    balanced_union = set().union(*(balanced[run] for run in common))
    result = {
        "status": "completed",
        "schema_version": "patient_balanced_regulatory_comparison.v1",
        "lambda_nominal": 0.001,
        "common_run_count": len(common),
        "pooled_edge_union": len(pooled_union),
        "patient_balanced_edge_union": len(balanced_union),
        "union_jaccard": _jaccard(pooled_union, balanced_union),
        "mean_sample_jaccard": sum(_jaccard(pooled[run], balanced[run]) for run in common) / len(common),
        "pooled_input": {"path": str(args.pooled), "sha256": _sha(args.pooled)},
        "patient_balanced_input": {"path": str(args.patient_balanced), "sha256": _sha(args.patient_balanced)},
        "response_blind": True,
        "claim_limit": "patient-weighting sensitivity only; no phenotype or causal interpretation",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "union_jaccard": result["union_jaccard"], "mean_sample_jaccard": result["mean_sample_jaccard"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
