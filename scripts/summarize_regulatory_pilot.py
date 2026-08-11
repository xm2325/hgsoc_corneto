** WARNING: connection is not using a post-quantum key exchange algorithm.
** This session may be vulnerable to "store now, decrypt later" attacks.
** The server may need to be upgraded. See https://openssh.com/pq.html
#!/usr/bin/env python3
"""Summarize a response-blind regulatory CORNETO receipt.

The output is descriptive only: sample status, edge/TF frequency and, when
requested, a run-accession join to the existing technical NMF rank-3 labels.
No phenotype values are read or inferred.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nmf-rank3", type=Path)
    args = parser.parse_args()
    pilot = json.loads(args.pilot.read_text(encoding="utf-8"))
    samples = pilot.get("samples")
    if not isinstance(samples, list) or not samples:
        raise ValueError("pilot receipt has no samples")
    status_counts = Counter(str(row.get("status", "missing")) for row in samples)
    edge_counts: Counter[tuple[str, str, int]] = Counter()
    tf_counts: Counter[str] = Counter()
    input_counts: Counter[str] = Counter()
    for row in samples:
        for edge in row.get("selected_edges", []):
            edge_counts[(str(edge.get("source")), str(edge.get("target")), int(edge.get("sign", 0)))] += 1
        for tf in row.get("outputs", {}):
            tf_counts[str(tf)] += 1
        for node in row.get("inputs", {}):
            input_counts[str(node)] += 1
    summary: dict[str, object] = {
        "schema_version": "regulatory_summary.v1",
        "response_blind": True,
        "pilot": str(args.pilot),
        "pilot_status": pilot.get("status"),
        "study": pilot.get("study"),
        "solver": pilot.get("solver"),
        "source_sha256": pilot.get("source_sha256", {}),
        "sample_count": len(samples),
        "status_counts": dict(sorted(status_counts.items())),
        "optimal_sample_count": sum(1 for row in samples if row.get("status") in {"optimal", "optimal_inaccurate"}),
        "edge_union_count": len(edge_counts),
        "edge_frequency": [
            {"source": source, "target": target, "sign": sign, "sample_count": count}
            for (source, target, sign), count in sorted(edge_counts.items(), key=lambda item: (-item[1], item[0]))
        ],
        "tf_frequency": [{"tf": tf, "sample_count": count} for tf, count in tf_counts.most_common()],
        "input_frequency": [{"node": node, "sample_count": count} for node, count in input_counts.most_common()],
        "nmf_integration": {"status": "not_requested"},
    }
    if args.nmf_rank3:
        rows = _read_tsv(args.nmf_rank3)
        by_run = {row.get("run_accession", ""): row for row in rows}
        groups: dict[str, list[dict[str, object]]] = defaultdict(list)
        for sample in samples:
            run = str(sample.get("run_accession", ""))
            nmf = by_run.get(run)
            if nmf is None:
                continue
            groups[nmf.get("independent_state", "unknown")].append({
                "run_accession": run,
                "regulatory_status": sample.get("status"),
                "selected_edge_count": len(sample.get("selected_edges", [])),
                "objective": sample.get("objective"),
                "assignment_stability": nmf.get("assignment_stability"),
            })
        summary["nmf_integration"] = {
            "status": "descriptive_run_join",
            "rank3_path": str(args.nmf_rank3),
            "matched_runs": sum(len(value) for value in groups.values()),
            "state_summary": {
                state: {
                    "sample_count": len(items),
                    "status_counts": dict(Counter(str(item["regulatory_status"]) for item in items)),
                    "mean_selected_edge_count": sum(float(item["selected_edge_count"]) for item in items) / len(items),
                    "runs": items,
                }
                for state, items in sorted(groups.items())
            },
            "claim_limit": "technical NMF state labels are not Barnes labels and are not causal/phenotypic states",
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
