#!/usr/bin/env python3
"""Consensus NMF for one primary-only cohort or the pooled 60-OCM cohort."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from hgsoc_corneto.nmf import run_consensus_nmf, select_top_mad


ALL_STUDIES = ("E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568")


def read_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:2] != ["gene_id", "gene_name"]:
            raise ValueError(f"unexpected matrix header: {path}")
        genes, values = [], []
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"wrong field count: {path}")
            genes.append(row[0])
            values.append([float(value) for value in row[2:]])
    return genes, header[2:], np.asarray(values, dtype=float)


def read_qc(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["run_accession"]: row for row in rows}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_matrix(path: Path, row_field: str, rows: list[str], columns: list[str], values: np.ndarray) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_field, *columns])
        for name, values_row in zip(rows, values, strict=True):
            writer.writerow([name, *(format(float(value), ".12g") for value in values_row)])


def load_primary(repo: Path, matrix_root: Path, studies: tuple[str, ...]):
    shared_genes: list[str] | None = None
    matrices, metadata = [], []
    for study in studies:
        genes, samples, matrix = read_matrix(matrix_root / study / "gene_log1p_tpm.tsv.gz")
        if shared_genes is None:
            shared_genes = genes
        elif genes != shared_genes:
            raise ValueError(f"gene identifiers/order differ for {study}")
        qc = read_qc(repo / "data/processed/rna" / study / "aggregation/sample_qc.tsv")
        if set(samples) != set(qc):
            raise ValueError(f"matrix/QC run mismatch for {study}")
        indices = [i for i, run in enumerate(samples) if qc[run]["primary_cohort_eligible"].lower() == "true"]
        matrices.append(matrix[:, indices])
        for index in indices:
            row = dict(qc[samples[index]])
            if row["sample_class"] != "tumour" or row["histotype_group"] != "HGSOC":
                raise ValueError(f"invalid primary sample: {samples[index]}")
            metadata.append(row)
    assert shared_genes is not None
    return shared_genes, metadata, np.concatenate(matrices, axis=1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--matrix-root", type=Path, required=True)
    parser.add_argument("--target", required=True, choices=(*ALL_STUDIES, "pooled60"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-genes", type=int, default=6000)
    parser.add_argument("--starts", type=int, default=100)
    parser.add_argument("--max-rank", type=int, default=10)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError(f"refusing to overwrite: {args.output_dir}")
    studies = ALL_STUDIES if args.target == "pooled60" else (args.target,)
    genes, metadata, matrix = load_primary(args.repo, args.matrix_root, studies)
    samples = [row["run_accession"] for row in metadata]
    expected = 60 if args.target == "pooled60" else {"E-MTAB-7223": 9, "E-MTAB-10801": 13, "E-MTAB-11000": 11, "E-MTAB-14568": 27}[args.target]
    if len(samples) != expected or matrix.shape[1] != expected:
        raise ValueError(f"{args.target}: expected {expected} primary samples, got {len(samples)}")
    selected_indices, mad = select_top_mad(matrix, genes, min(args.top_genes, len(genes)))
    selected_genes = [genes[index] for index in selected_indices]
    selected = matrix[selected_indices, :]
    ranks = list(range(2, min(args.max_rank, len(samples) - 1) + 1))
    seed_base = 6000000 if args.target == "pooled60" else sum((i + 1) * ord(c) for i, c in enumerate(args.target)) * 100
    args.output_dir.mkdir(parents=True)
    write_tsv(args.output_dir / "selected_genes.tsv", [
        {"selection_rank": rank, "gene_id": genes[index], "mad": format(float(mad[index]), ".12g")}
        for rank, index in enumerate(selected_indices, start=1)
    ])
    summaries = []
    for rank in ranks:
        result = run_consensus_nmf(
            selected,
            tuple(samples),
            rank=rank,
            runs=args.starts,
            seed_base=seed_base + rank * args.starts,
            max_iter=2000,
            tolerance=1e-4,
        )
        rank_dir = args.output_dir / f"rank_{rank}"
        rank_dir.mkdir()
        write_matrix(rank_dir / "best_w_gene_loadings.tsv.gz", "gene_id", selected_genes, list(result.state_names), result.best_w)
        write_matrix(rank_dir / "best_h_sample_loadings.tsv.gz", "state", list(result.state_names), samples, result.best_h)
        write_matrix(rank_dir / "consensus.tsv.gz", "run_accession", samples, samples, result.consensus)
        assignments = []
        for index, row in enumerate(metadata):
            state = int(result.labels[index])
            assignments.append({
                "study_accession": row["study_accession"],
                "run_accession": row["run_accession"],
                "canonical_ocm_id": row["canonical_ocm_id"],
                "patient_id": row["patient_id"],
                "state": result.state_names[state],
                "assignment_stability": format(float(result.assignment_stability[index]), ".12g"),
                "consensus_silhouette": format(float(result.silhouette_by_sample[index]), ".12g"),
            })
        write_tsv(rank_dir / "assignments.tsv", assignments)
        best = min(result.runs, key=lambda run: run.reconstruction_error)
        fits = np.asarray([run.fit for run in result.runs])
        summaries.append({
            "rank": rank,
            "sample_count": len(samples),
            "gene_count": len(selected_genes),
            "nmf_runs": len(result.runs),
            "best_fit": format(float(best.fit), ".12g"),
            "median_fit": format(float(np.median(fits)), ".12g"),
            "cophenetic_correlation": format(float(result.cophenetic_correlation), ".12g"),
            "dispersion": format(float(result.dispersion), ".12g"),
            "consensus_sharpness": format(float(result.consensus_sharpness), ".12g"),
            "average_silhouette": format(float(result.average_silhouette), ".12g"),
            "cluster_sizes": ";".join(str(int(np.sum(result.labels == state))) for state in range(rank)),
            "converged_runs": sum(run.converged for run in result.runs),
        })
    write_tsv(args.output_dir / "rank_summary.tsv", summaries)
    receipt = {
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "target": args.target,
        "studies": list(studies),
        "sample_policy": "primary_cohort_eligible=true; tumour; HGSOC",
        "sample_count": len(samples),
        "patient_count": len({row["patient_id"] for row in metadata}),
        "gene_count": len(selected_genes),
        "ranks": ranks,
        "starts_per_rank": args.starts,
        "input_value": "log1p_tpm",
        "claim_limit": "response-blind technical states; pooled states must be checked for study association",
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    (args.output_dir / "nmf_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
