#!/usr/bin/env python3
"""Compare narrow and richer true-joint PKN receipts at fixed lambda."""

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
    result = {}
    for row in rows:
        run = row.get("run_accession")
        if not isinstance(run, str) or run in result:
            raise ValueError(f"{path} has invalid/duplicate run")
        edges = {(edge["source"], edge["target"], int(edge["sign"])) for edge in row.get("selected_edges", [])}
        result[run] = edges
    return result


def _jaccard(left: set, right: set) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def compare(narrow_root: Path, richer_root: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise ValueError(f"refusing to overwrite {output}")
    groups = ["pooled", "E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568"]
    rows = []
    for group in groups:
        narrow = narrow_root / group / "l0p001.json"
        richer = richer_root / group / "l0p001.json"
        left, right = _read(narrow), _read(richer)
        common = sorted(set(left) & set(right))
        if set(left) != set(right):
            raise ValueError(f"run set mismatch for {group}")
        left_union = set().union(*(left[run] for run in common))
        right_union = set().union(*(right[run] for run in common))
        rows.append(
            {
                "group": group,
                "sample_count": len(common),
                "narrow_edge_union": len(left_union),
                "richer_edge_union": len(right_union),
                "union_jaccard": _jaccard(left_union, right_union),
                "mean_sample_jaccard": sum(_jaccard(left[run], right[run]) for run in common) / len(common),
                "narrow_sha256": _sha(narrow),
                "richer_sha256": _sha(richer),
            }
        )
    result = {
        "status": "completed",
        "schema_version": "regulatory_richer_pkn_comparison.v1",
        "lambda_nominal": 0.001,
        "narrow_policy": "v1 max_inputs=3 max_outputs=6 max_depth=3",
        "richer_policy": "richer_v1 max_inputs=5 max_outputs=10 max_depth=4",
        "groups": rows,
        "response_blind": True,
        "claim_limit": "PKN/graph-policy sensitivity only; no phenotype or causal interpretation",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--narrow-root", type=Path, required=True)
    parser.add_argument("--richer-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.narrow_root, args.richer_root, args.output)
    print(json.dumps({"status": result["status"], "groups": len(result["groups"])}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
