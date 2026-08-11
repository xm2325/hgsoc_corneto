#!/usr/bin/env python3
"""Run the frozen independent E-MTAB-14568 consensus-NMF benchmark."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import platform
import subprocess
import tempfile
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import scipy
import sklearn
import yaml

from hgsoc_corneto.io import (
    deterministic_gzip_text_writer,
    read_tsv,
    sha256,
    write_json,
    write_tsv,
)
from hgsoc_corneto.nmf import ConsensusNmfResult, run_consensus_nmf, select_top_mad


def _git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, text=True, stdout=subprocess.PIPE
    )
    return result.stdout.strip()


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open()


def _read_matrix(path: Path) -> tuple[tuple[str, ...], tuple[str, ...], np.ndarray]:
    gene_ids: list[str] = []
    rows: list[list[float]] = []
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"Empty expression matrix: {path}") from error
        if len(header) < 4 or header[:2] != ["gene_id", "gene_name"]:
            raise ValueError(f"Unexpected expression matrix header: {path}")
        sample_ids = tuple(header[2:])
        if len(set(sample_ids)) != len(sample_ids):
            raise ValueError("Expression matrix sample identifiers must be unique")
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Wrong field count at {path}:{line_number}")
            gene_ids.append(row[0])
            rows.append([float(value) for value in row[2:]])
    if len(set(gene_ids)) != len(gene_ids):
        raise ValueError("Expression matrix gene identifiers must be unique")
    matrix = np.asarray(rows, dtype=float)
    if matrix.shape != (len(gene_ids), len(sample_ids)):
        raise ValueError("Expression matrix shape mismatch")
    return tuple(gene_ids), sample_ids, matrix


def _write_numeric_matrix(
    path: Path,
    row_name: str,
    row_ids: tuple[str, ...],
    column_ids: tuple[str, ...],
    values: np.ndarray,
) -> None:
    with deterministic_gzip_text_writer(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_name, *column_ids])
        for row_id, row in zip(row_ids, values, strict=True):
            writer.writerow([row_id, *(format(float(value), ".12g") for value in row)])


def _rank_outputs(
    output_dir: Path,
    result: ConsensusNmfResult,
    sample_ids: tuple[str, ...],
    selected_gene_ids: tuple[str, ...],
) -> None:
    rank_dir = output_dir / f"rank_{result.rank}"
    rank_dir.mkdir()
    _write_numeric_matrix(
        rank_dir / "consensus.tsv.gz",
        "run_accession",
        sample_ids,
        sample_ids,
        result.consensus,
    )
    _write_numeric_matrix(
        rank_dir / "best_w_gene_loadings.tsv.gz",
        "gene_id",
        selected_gene_ids,
        result.state_names,
        result.best_w,
    )
    _write_numeric_matrix(
        rank_dir / "best_h_sample_loadings.tsv.gz",
        "state",
        result.state_names,
        sample_ids,
        result.best_h,
    )
    write_tsv(
        rank_dir / "run_metrics.tsv",
        [
            {
                "run_index": run.run_index,
                "seed": run.seed,
                "reconstruction_error": format(run.reconstruction_error, ".12g"),
                "fit": format(run.fit, ".12g"),
                "iterations": run.iterations,
                "converged": run.converged,
            }
            for run in result.runs
        ],
    )
    with deterministic_gzip_text_writer(rank_dir / "aligned_run_assignments.tsv.gz") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["run_accession", *(f"nmf_run_{run.run_index}" for run in result.runs)])
        for sample_index, sample_id in enumerate(sample_ids):
            writer.writerow(
                [
                    sample_id,
                    *(
                        result.state_names[int(run.labels[sample_index])]
                        for run in result.runs
                    ),
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--sample-qc", type=Path, required=True)
    parser.add_argument("--artifact-output-dir", type=Path, required=True)
    parser.add_argument("--report-output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.artifact_output_dir.exists() or args.report_output_dir.exists():
        raise ValueError("Existing NMF output requires inspection; refusing to overwrite")
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    analysis = config["independent_etab_14568"]
    gene_ids, sample_ids, matrix = _read_matrix(args.matrix)
    qc_rows = read_tsv(args.sample_qc)
    qc_by_run = {row["run_accession"]: row for row in qc_rows}
    if set(sample_ids) != set(qc_by_run) or len(qc_rows) != len(qc_by_run):
        raise ValueError("Expression matrix and sample QC disagree")
    if len(sample_ids) != 33 or analysis["sample_policy"] != "all_33_public_runs":
        raise ValueError("Frozen E-MTAB-14568 benchmark requires all 33 public runs")

    selected_indices, mad = select_top_mad(matrix, gene_ids, int(analysis["top_genes"]))
    selected_gene_ids = tuple(gene_ids[index] for index in selected_indices)
    selected_matrix = matrix[selected_indices, :]
    implementation = analysis["implementation"]

    args.artifact_output_dir.parent.mkdir(parents=True, exist_ok=True)
    args.report_output_dir.parent.mkdir(parents=True, exist_ok=True)
    artifact_stage = Path(
        tempfile.mkdtemp(prefix=".nmf-artifacts.", dir=args.artifact_output_dir.parent)
    )
    report_stage = Path(
        tempfile.mkdtemp(prefix=".nmf-report.", dir=args.report_output_dir.parent)
    )
    try:
        write_tsv(
            artifact_stage / "selected_genes.tsv",
            [
                {
                    "selection_rank": rank,
                    "gene_id": gene_ids[index],
                    "mad": format(float(mad[index]), ".12g"),
                }
                for rank, index in enumerate(selected_indices, start=1)
            ],
        )
        results: dict[int, ConsensusNmfResult] = {}
        summary_rows = []
        for rank in (int(value) for value in analysis["ranks"]):
            result = run_consensus_nmf(
                selected_matrix,
                sample_ids,
                rank=rank,
                runs=int(analysis["random_initializations_per_rank"]),
                seed_base=int(analysis["seed_base"]),
                max_iter=int(implementation["maximum_iterations"]),
                tolerance=float(implementation["tolerance"]),
            )
            results[rank] = result
            _rank_outputs(artifact_stage, result, sample_ids, selected_gene_ids)
            best_run = min(result.runs, key=lambda run: run.reconstruction_error)
            fits = np.asarray([run.fit for run in result.runs])
            cluster_sizes = [int(np.sum(result.labels == state)) for state in range(rank)]
            summary_rows.append(
                {
                    "rank": rank,
                    "sample_count": len(sample_ids),
                    "gene_count": len(selected_gene_ids),
                    "nmf_runs": len(result.runs),
                    "best_seed": best_run.seed,
                    "best_reconstruction_error": format(best_run.reconstruction_error, ".12g"),
                    "best_fit": format(best_run.fit, ".12g"),
                    "median_fit": format(float(np.median(fits)), ".12g"),
                    "cophenetic_correlation": format(result.cophenetic_correlation, ".12g"),
                    "dispersion": format(result.dispersion, ".12g"),
                    "consensus_sharpness": format(result.consensus_sharpness, ".12g"),
                    "average_silhouette": format(result.average_silhouette, ".12g"),
                    "cluster_sizes": ";".join(str(value) for value in cluster_sizes),
                    "converged_runs": sum(run.converged for run in result.runs),
                }
            )
        write_tsv(report_stage / "rank_summary.tsv", summary_rows)

        primary_rank = int(analysis["primary_rank"])
        primary = results[primary_rank]
        normalized_h = primary.best_h / np.maximum(primary.best_h.sum(axis=0), 1e-15)
        assignments = []
        for sample_index, sample_id in enumerate(sample_ids):
            qc = qc_by_run[sample_id]
            state_index = int(primary.labels[sample_index])
            assignments.append(
                {
                    "study_accession": analysis["study_accession"],
                    "run_accession": sample_id,
                    "canonical_ocm_id": qc["canonical_ocm_id"],
                    "patient_id": qc["patient_id"],
                    "histotype_group": qc["histotype_group"],
                    "primary_cohort_eligible": qc["primary_cohort_eligible"],
                    "independent_state": primary.state_names[state_index],
                    "assignment_stability": format(
                        float(primary.assignment_stability[sample_index]), ".12g"
                    ),
                    "consensus_silhouette": format(
                        float(primary.silhouette_by_sample[sample_index]), ".12g"
                    ),
                    "best_solution_state_loading_fraction": format(
                        float(normalized_h[state_index, sample_index]), ".12g"
                    ),
                }
            )
        write_tsv(report_stage / "rank_3_assignments.tsv", assignments)

        artifact_files = sorted(path for path in artifact_stage.rglob("*") if path.is_file())
        receipt = {
            "status": "completed",
            "completed_at_utc": datetime.now(UTC).isoformat(),
            "repo_commit": _git_commit(),
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "scikit_learn": sklearn.__version__,
            },
            "config": config,
            "input": {
                "matrix_path": str(args.matrix),
                "matrix_sha256": sha256(args.matrix),
                "sample_qc_path": str(args.sample_qc),
                "sample_qc_sha256": sha256(args.sample_qc),
                "samples": len(sample_ids),
                "genes_before_mad_selection": len(gene_ids),
                "genes_after_mad_selection": len(selected_gene_ids),
            },
            "artifacts": [
                {
                    "path": str(args.artifact_output_dir / path.relative_to(artifact_stage)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in artifact_files
            ],
            "reports": {
                "rank_summary": str(args.report_output_dir / "rank_summary.tsv"),
                "rank_3_assignments": str(args.report_output_dir / "rank_3_assignments.tsv"),
            },
            "claim_limit": config["barnes_method_boundary"]["claim_limit"],
        }
        write_json(report_stage / "nmf_receipt.json", receipt)
        os.replace(artifact_stage, args.artifact_output_dir)
        os.replace(report_stage, args.report_output_dir)
    except Exception as error:
        raise RuntimeError(
            f"NMF staging retained: artifacts={artifact_stage}; report={report_stage}"
        ) from error
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
