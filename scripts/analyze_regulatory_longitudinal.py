** WARNING: connection is not using a post-quantum key exchange algorithm.
** This session may be vulnerable to "store now, decrypt later" attacks.
** The server may need to be upgraded. See https://openssh.com/pq.html
#!/usr/bin/env python3
"""Response-blind within-family regulatory edge-change summary."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def _split(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pilot = json.loads(args.pilot.read_text(encoding="utf-8"))
    sample_by_run = {row.get("run_accession"): row for row in pilot.get("samples", [])}
    with args.families.open(encoding="utf-8", newline="") as handle:
        families = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for family in families:
        runs = _split(family.get("primary_run_accessions", ""))
        available = [sample_by_run[run] for run in runs if run in sample_by_run]
        if len(available) < 2:
            continue
        baseline = available[0]
        baseline_edges = {(e.get("source"), e.get("target"), int(e.get("sign", 0))) for e in baseline.get("selected_edges", [])}
        for current in available[1:]:
            current_edges = {(e.get("source"), e.get("target"), int(e.get("sign", 0))) for e in current.get("selected_edges", [])}
            union = baseline_edges | current_edges
            rows.append({
                "family_id": family.get("family_id", family.get("patient_id", "")),
                "patient_id": family.get("patient_id", ""),
                "analysis_role": family.get("analysis_role", ""),
                "pair_status": family.get("pair_status", ""),
                "baseline_run": baseline.get("run_accession"),
                "current_run": current.get("run_accession"),
                "baseline_status": baseline.get("status"),
                "current_status": current.get("status"),
                "baseline_edge_count": len(baseline_edges),
                "current_edge_count": len(current_edges),
                "gained_edge_count": len(current_edges - baseline_edges),
                "lost_edge_count": len(baseline_edges - current_edges),
                "edge_jaccard": len(baseline_edges & current_edges) / len(union) if union else 1.0,
            })
    result = {
        "schema_version": "regulatory_longitudinal.v1",
        "response_blind": True,
        "pilot": str(args.pilot),
        "families": str(args.families),
        "pair_count": len(rows),
        "status_counts": dict(Counter(
            f"{row['baseline_status']}__{row['current_status']}" for row in rows
        )),
        "pairs": rows,
        "claim_limit": "within-family descriptive rewiring only; no treatment causality or phenotype association",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
