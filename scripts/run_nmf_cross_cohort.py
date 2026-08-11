#!/usr/bin/env python3
"""Run the response-blind consensus-NMF benchmark on one non-14568 cohort.

This deliberately keeps the frozen 14568 settings (log1p TPM, top 6000 MAD
genes, ranks 2--8, 100 random starts) while allowing the cohort sample count
to vary.  It writes a compact rank summary and rank-3 run assignments; the
labels are technical cohort-local states and are not Barnes subtype labels.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from hgsoc_corneto.nmf import run_consensus_nmf, select_top_mad


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with _open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:2] != ["gene_id", "gene_name"]:
            raise ValueError(f"unexpected matrix header: {path}")
        samples = header[2:]
        genes: list[str] = []
        values: list[list[float]] = []
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"wrong field count in {path}")
            genes.append(row[0])
            values.append([float(value) for value in row[2:]])
    matrix = np.asarray(values, dtype=float)
    if matrix.shape != (len(genes), len(samples)):
        raise ValueError("matrix shape mismatch")
    return genes, samples, matrix


def _qc(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["run_accession"]: row for row in rows}


def _write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--sample-qc", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-genes", type=int, default=6000)
    parser.add_argument("--starts", type=int, default=100)
    parser.add_argument("--max-rank", type=int, default=8)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError(f"refusing to overwrite {args.output_dir}")
    genes, samples, matrix = _matrix(args.matrix)
    qc = _qc(args.sample_qc)
    if set(samples) != set(qc):
        raise ValueError("matrix and sample QC run sets disagree")
    if len(samples) < 4:
        raise ValueError("too few samples for cross-cohort NMF")
    selected_indices, mad = select_top_mad(matrix, genes, min(args.top_genes, len(genes)))
    selected = matrix[selected_indices, :]
    ranks = list(range(2, min(args.max_rank, len(samples) - 1) + 1))
    seed_base = sum((index + 1) * ord(char) for index, char in enumerate(args.study)) * 100
    summaries: list[dict[str, object]] = []
    rank3 = None
    for rank in ranks:
        result = run_consensus_nmf(
            selected,
            tuple(samples),
            rank=rank,
            runs=args.starts,
            seed_base=seed_base + rank,
            max_iter=2000,
            tolerance=1e-4,
        )
        best = min(result.runs, key=lambda run: run.reconstruction_error)
        fits = np.asarray([run.fit for run in result.runs])
        summaries.append({
            "rank": rank,
            "sample_count": len(samples),
            "gene_count": len(selected_indices),
            "nmf_runs": len(result.runs),
            "best_fit": float(best.fit),
            "median_fit": float(np.median(fits)),
            "cophenetic_correlation": float(result.cophenetic_correlation),
            "average_silhouette": float(result.average_silhouette),
            "cluster_sizes": [int(np.sum(result.labels == state)) for state in range(rank)],
        })
        if rank == 3:
            rank3 = result
    if rank3 is None:
        raise ValueError("rank 3 is unavailable for this cohort")
    args.output_dir.mkdir(parents=True)
    assignments = []
    for index, run in enumerate(samples):
        row = qc[run]
        state_index = int(rank3.labels[index])
        assignments.append({
            "study_accession": args.study,
            "run_accession": run,
            "canonical_ocm_id": row.get("canonical_ocm_id", ""),
            "patient_id": row.get("patient_id", ""),
            "primary_cohort_eligible": row.get("primary_cohort_eligible", ""),
            "independent_state": rank3.state_names[state_index],
            "assignment_stability": float(rank3.assignment_stability[index]),
            "consensus_silhouette": float(rank3.silhouette_by_sample[index]),
        })
    _write_tsv(args.output_dir / "rank_3_assignments.tsv", assignments, list(assignments[0]))
    _write_tsv(args.output_dir / "rank_summary.tsv", summaries, list(summaries[0]))
    receipt = {
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "study_accession": args.study,
        "sample_count": len(samples),
        "gene_count": len(selected_indices),
        "input_value": "log1p_tpm",
        "top_genes": min(args.top_genes, len(genes)),
        "ranks": ranks,
        "random_initializations_per_rank": args.starts,
        "rank3_assignment_count": len(assignments),
        "claim_limit": "technical cohort-local NMF states; no Barnes subtype or drug-response interpretation",
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    (args.output_dir / "nmf_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "study": args.study, "samples": len(samples), "output": str(args.output_dir)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
