#!/usr/bin/env python3
"""Fail-closed aggregation of response-blind regulatory solution ensembles."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _median(values: list[float | int | None]) -> float | None:
    clean = [float(value) for value in values if value is not None]
    return None if not clean else statistics.median(clean)


def summarize(paths: list[Path], expected_samples: int, study: str) -> dict[str, Any]:
    if len(paths) < expected_samples:
        raise ValueError(f"found {len(paths)} receipts, expected at least {expected_samples}")
    candidates: dict[str, list[dict[str, Any]]] = {}
    method_fingerprints = set()
    source_fingerprints = set()
    receipt_files = []
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema_version") != "regulatory_alternative_optima.v1":
            raise ValueError(f"{path}: unexpected schema")
        if payload.get("study") != study:
            raise ValueError(f"{path}: study mismatch")
        rows = payload.get("samples")
        if not isinstance(rows, list) or len(rows) != 1:
            raise ValueError(f"{path}: expected exactly one sample")
        sample = dict(rows[0])
        run = str(sample.get("run_accession", ""))
        if not run:
            raise ValueError(f"{path}: missing run_accession")
        # Older raw receipts may predate the explicit fail-closed status.  Keep
        # the raw status but classify an all-zero incumbent as blocked here.
        raw_status = str(sample.get("status", "unknown"))
        status = raw_status
        if raw_status == "completed" and sample.get("incumbent_edge_count") == 0:
            status = "blocked_no_selected_edges"
        candidates.setdefault(run, []).append(
            {
                "run_accession": run,
                "status": status,
                "raw_status": raw_status,
                "selected_receipt": str(path),
                "solution_count": sample.get("solution_count"),
                "accepted_alternative_count": sample.get("accepted_alternative_count"),
                "incumbent_edge_count": sample.get("incumbent_edge_count"),
                "edge_union_count": sample.get("edge_union_count"),
                "edge_intersection_count": sample.get("edge_intersection_count"),
                "core_edge_count": sample.get("core_edge_count"),
                "mean_pairwise_jaccard": sample.get("mean_pairwise_jaccard"),
                "min_pairwise_jaccard": sample.get("min_pairwise_jaccard"),
                "mean_incumbent_jaccard": sample.get("mean_incumbent_jaccard"),
                "mean_edge_entropy_bits": sample.get("mean_edge_entropy_bits"),
            }
        )
        method_fingerprints.add(json.dumps(payload.get("method"), sort_keys=True))
        source_fingerprints.add(json.dumps(payload.get("source_sha256"), sort_keys=True))
        receipt_files.append({"path": str(path), "sha256": _sha256(path)})
    if len(method_fingerprints) != 1:
        raise ValueError("receipts disagree on method parameters")
    if len(source_fingerprints) != 1:
        raise ValueError("receipts disagree on source hashes")
    if len(candidates) != expected_samples:
        raise ValueError(
            f"receipts represent {len(candidates)} unique samples, expected {expected_samples}"
        )

    def _candidate_rank(row: dict[str, Any]) -> tuple[int, str]:
        status = str(row["status"])
        if status == "completed":
            rank = 3
        elif status.startswith("blocked"):
            rank = 2
        else:
            rank = 1
        return rank, str(row["selected_receipt"])

    samples = [max(rows, key=_candidate_rank) for rows in candidates.values()]
    samples.sort(key=lambda row: row["run_accession"])
    usable = [row for row in samples if row["status"] == "completed"]
    status_counts = Counter(row["status"] for row in samples)
    return {
        "schema_version": "regulatory_alternative_optima_summary.v1",
        "status": "completed" if usable else "blocked",
        "created_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "response_blind": True,
        "phenotype_inputs": [],
        "study": study,
        "method": json.loads(next(iter(method_fingerprints))),
        "source_sha256": json.loads(next(iter(source_fingerprints))),
        "counts": {
            "expected_samples": expected_samples,
            "input_receipts": len(paths),
            "receipt_samples": len(samples),
            "superseded_receipts": len(paths) - len(samples),
            "usable_nonempty_samples": len(usable),
            "status_counts": dict(sorted(status_counts.items())),
        },
        "ensemble_metrics_nonempty_samples": {
            "median_accepted_alternative_count": _median(
                [row["accepted_alternative_count"] for row in usable]
            ),
            "median_incumbent_edge_count": _median(
                [row["incumbent_edge_count"] for row in usable]
            ),
            "median_edge_union_count": _median([row["edge_union_count"] for row in usable]),
            "median_core_edge_count": _median([row["core_edge_count"] for row in usable]),
            "median_mean_pairwise_jaccard": _median(
                [row["mean_pairwise_jaccard"] for row in usable]
            ),
            "median_min_pairwise_jaccard": _median(
                [row["min_pairwise_jaccard"] for row in usable]
            ),
            "median_mean_edge_entropy_bits": _median(
                [row["mean_edge_entropy_bits"] for row in usable]
            ),
        },
        "samples": samples,
        "receipts": receipt_files,
        "interpretation_guardrail": (
            "Near-optimal CARNIVAL edge stability is response-blind and does not establish "
            "Taxol association, drug resistance, or causality. All-zero incumbents are blocked."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, action="append", required=True)
    parser.add_argument("--pattern", default="sample_*.json")
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = sorted(path for directory in args.input_dir for path in directory.glob(args.pattern))
    result = summarize(paths, args.expected_samples, args.study)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
