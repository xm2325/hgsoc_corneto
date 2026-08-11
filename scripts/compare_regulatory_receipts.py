#!/usr/bin/env python3
"""Compare response-blind regulatory receipts across cohorts."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="append", required=True, metavar="STUDY=JSON")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cohorts = {}
    for item in args.pilot:
        if "=" not in item:
            raise ValueError("--pilot must be STUDY=JSON")
        study, path = item.split("=", 1)
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data.get("status") != "completed":
            raise ValueError(f"{study}: receipt status is {data.get('status')!r}")
        edges = {(e.get("source"), e.get("target"), int(e.get("sign", 0)))
                 for sample in data.get("samples", [])
                 for e in sample.get("selected_edges", [])}
        tfs = {tf for sample in data.get("samples", []) for tf in sample.get("outputs", {})}
        cohorts[study] = {"path": path, "samples": len(data.get("samples", [])),
                          "status_counts": dict(Counter(s.get("status") for s in data.get("samples", []))),
                          "edges": edges, "tfs": tfs}
    pairwise = []
    studies = sorted(cohorts)
    for index, left in enumerate(studies):
        for right in studies[index + 1:]:
            a, b = cohorts[left]["edges"], cohorts[right]["edges"]
            union = a | b
            pairwise.append({"left": left, "right": right, "intersection": len(a & b),
                             "union": len(union), "jaccard": (len(a & b) / len(union)) if union else 1.0})
    edge_frequency = Counter(edge for value in cohorts.values() for edge in value["edges"])
    result = {
        "schema_version": "regulatory_cross_cohort.v1", "response_blind": True,
        "cohorts": {study: {k: v for k, v in value.items() if k not in {"edges", "tfs"}}
                    | {"edge_count": len(value["edges"]), "tf_count": len(value["tfs"])}
                    for study, value in sorted(cohorts.items())},
        "pairwise_edge_jaccard": pairwise,
        "edge_cohort_frequency": [{"source": edge[0], "target": edge[1], "sign": edge[2], "cohort_count": count}
                                  for edge, count in edge_frequency.most_common()],
        "claim_limit": "descriptive response-blind overlap; no drug-response or causal interpretation",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
