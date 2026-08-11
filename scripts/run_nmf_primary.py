#!/usr/bin/env python3
"""Consensus NMF for primary tumour OCMs in one cohort or a pooled cohort.

The runner deliberately writes to a new output namespace.  It never modifies
the historical all-run E-MTAB-14568 benchmark.  Every input cohort is filtered
to ``primary_cohort_eligible=true`` using its aggregation QC table, and pooled
analyses retain the originating study in every sample-level output.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import platform
import subprocess
import tempfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import scipy
import sklearn
from scipy.stats import chi2_contingency

from hgsoc_corneto.io import deterministic_gzip_text_writer, sha256, write_json, write_tsv
from hgsoc_corneto.nmf import ConsensusNmfResult, run_consensus_nmf, select_top_mad


def _git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, text=True, stdout=subprocess.PIPE
    )
    return result.stdout.strip()


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _read_matrix(path: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], np.ndarray]:
    gene_ids: list[str] = []
    gene_names: list[str] = []
    rows: list[list[float]] = []
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 4 or header[:2] != ["gene_id", "gene_name"]:
            raise ValueError(f"Unexpected expression matrix header: {path}")
        sample_ids = tuple(header[2:])
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError(f"Duplicate sample identifiers in {path}")
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Wrong field count at {path}:{line_number}")
            gene_ids.append(row[0])
            gene_names.append(row[1])
            rows.append([float(value) for value in row[2:]])
    if len(gene_ids) != len(set(gene_ids)):
        raise ValueError(f"Duplicate gene identifiers in {path}")
    values = np.asarray(rows, dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"NMF matrix must be finite and non-negative: {path}")
    return tuple(gene_ids), tuple(gene_names), sample_ids, values


def _read_qc(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "study_accession",
        "run_accession",
        "canonical_ocm_id",
        "patient_id",
        "sample_class",
        "histotype_group",
        "primary_cohort_eligible",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Missing required sample-QC columns: {path}")
    result = {row["run_accession"]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate run accessions in {path}")
    return result


def _parse_spec(value: str) -> tuple[str, Path, Path]:
    fields = value.split("::")
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("cohort spec must be STUDY::MATRIX::SAMPLE_QC")
    return fields[0], Path(fields[1]), Path(fields[2])


def load_primary_cohorts(
    specs: list[tuple[str, Path, Path]],
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], np.ndarray, list[dict[str, str]], list[dict[str, object]]]:
    """Load and concatenate primary tumour samples, preserving study metadata."""

    if not specs:
        raise ValueError("At least one cohort is required")
    reference_gene_ids: tuple[str, ...] | None = None
    reference_gene_names: tuple[str, ...] | None = None
    all_samples: list[str] = []
    all_values: list[np.ndarray] = []
    all_metadata: list[dict[str, str]] = []
    inputs: list[dict[str, object]] = []
    seen_studies: set[str] = set()
    for study, matrix_path, qc_path in specs:
        if study in seen_studies:
            raise ValueError(f"Duplicate cohort spec: {study}")
        seen_studies.add(study)
        gene_ids, gene_names, matrix_samples, matrix = _read_matrix(matrix_path)
        qc = _read_qc(qc_path)
        if set(matrix_samples) != set(qc):
            raise ValueError(f"Matrix/QC sample mismatch for {study}")
        if any(row["study_accession"] != study for row in qc.values()):
            raise ValueError(f"QC study accession mismatch for {study}")
        if reference_gene_ids is None:
            reference_gene_ids = gene_ids
            reference_gene_names = gene_names
        elif gene_ids != reference_gene_ids or gene_names != reference_gene_names:
            raise ValueError(f"Gene reference/order differs for {study}; pooled NMF is unsafe")
        keep_indices = [
            index
            for index, run in enumerate(matrix_samples)
            if qc[run]["primary_cohort_eligible"].strip().lower() == "true"
        ]
        if not keep_indices:
            raise ValueError(f"No primary-cohort samples for {study}")
        kept_samples = [matrix_samples[index] for index in keep_indices]
        if any(qc[run]["sample_class"] != "tumour" or qc[run]["histotype_group"] != "HGSOC" for run in kept_samples):
            raise ValueError(f"Primary policy includes a non-HGSOC/non-tumour sample for {study}")
        all_samples.extend(kept_samples)
        all_values.append(matrix[:, keep_indices])
        all_metadata.extend(qc[run] for run in kept_samples)
        inputs.append(
            {
                "study_accession": study,
                "matrix_path": str(matrix_path),
                "matrix_sha256": sha256(matrix_path),
                "sample_qc_path": str(qc_path),
                "sample_qc_sha256": sha256(qc_path),
                "matrix_samples": len(matrix_samples),
                "primary_samples": len(kept_samples),
            }
        )
    if len(all_samples) != len(set(all_samples)):
        raise ValueError("Run accessions are not unique across pooled cohorts")
    assert reference_gene_ids is not None and reference_gene_names is not None
    return (
        reference_gene_ids,
        reference_gene_names,
        tuple(all_samples),
        np.concatenate(all_values, axis=1),
        all_metadata,
        inputs,
    )


def _write_numeric_matrix(
    path: Path,
    row_label: str,
    row_ids: tuple[str, ...],
    column_ids: tuple[str, ...],
    values: np.ndarray,
) -> None:
    with deterministic_gzip_text_writer(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_label, *column_ids])
        for row_id, row in zip(row_ids, values, strict=True):
            writer.writerow([row_id, *(format(float(value), ".12g") for value in row)])


def _study_association(labels: np.ndarray, studies: tuple[str, ...], rank: int) -> dict[str, object]:
    study_names = tuple(sorted(set(studies)))
    table = np.asarray(
        [[np.sum((np.asarray(studies) == study) & (labels == state)) for state in range(rank)] for study in study_names],
        dtype=int,
    )
    if len(study_names) == 1:
        return {"chi_square": 0.0, "p_value": 1.0, "cramers_v": 0.0, "table": table.tolist()}
    chi_square, p_value, _, _ = chi2_contingency(table)
    denominator = len(labels) * min(table.shape[0] - 1, table.shape[1] - 1)
    cramers_v = float(np.sqrt(chi_square / denominator)) if denominator > 0 else 0.0
    return {
        "chi_square": float(chi_square),
        "p_value": float(p_value),
        "cramers_v": cramers_v,
        "table": table.tolist(),
        "study_order": list(study_names),
    }


def _rank_outputs(
    output_dir: Path,
    result: ConsensusNmfResult,
    sample_ids: tuple[str, ...],
    metadata: list[dict[str, str]],
    selected_gene_ids: tuple[str, ...],
) -> None:
    rank_dir = output_dir / f"rank_{result.rank}"
    rank_dir.mkdir()
    _write_numeric_matrix(rank_dir / "consensus.tsv.gz", "run_accession", sample_ids, sample_ids, result.consensus)
    _write_numeric_matrix(
        rank_dir / "gene_loadings.tsv.gz", "gene_id", selected_gene_ids, result.state_names, result.best_w
    )
    _write_numeric_matrix(
        rank_dir / "sample_loadings.tsv.gz", "state", result.state_names, sample_ids, result.best_h
    )
    normalized_h = result.best_h / np.maximum(result.best_h.sum(axis=0), 1e-15)
    assignments: list[dict[str, object]] = []
    for sample_index, (sample_id, row) in enumerate(zip(sample_ids, metadata, strict=True)):
        state = int(result.labels[sample_index])
        assignments.append(
            {
                "study_accession": row["study_accession"],
                "run_accession": sample_id,
                "canonical_ocm_id": row["canonical_ocm_id"],
                "patient_id": row["patient_id"],
                "primary_cohort_eligible": row["primary_cohort_eligible"],
                "technical_state": result.state_names[state],
                "assignment_stability": format(float(result.assignment_stability[sample_index]), ".12g"),
                "consensus_silhouette": format(float(result.silhouette_by_sample[sample_index]), ".12g"),
                "best_solution_state_loading_fraction": format(float(normalized_h[state, sample_index]), ".12g"),
            }
        )
    write_tsv(rank_dir / "assignments.tsv", assignments)
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--cohort", action="append", type=_parse_spec, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-genes", type=int, default=6000)
    parser.add_argument("--ranks", default="2,3,4,5,6,7,8,9,10")
    parser.add_argument("--starts", type=int, default=100)
    parser.add_argument("--seed-base", type=int, default=608110)
    parser.add_argument("--max-iter", type=int, default=2000)
    parser.add_argument("--tolerance", type=float, default=1e-4)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError(f"Refusing to overwrite existing NMF output: {args.output_dir}")

    gene_ids, gene_names, samples, matrix, metadata, inputs = load_primary_cohorts(args.cohort)
    if len(samples) < 5:
        raise ValueError("At least five primary samples are required")
    requested_ranks = tuple(dict.fromkeys(int(value) for value in args.ranks.split(",")))
    ranks = tuple(rank for rank in requested_ranks if 2 <= rank <= len(samples) - 2)
    if 3 not in ranks:
        raise ValueError("Rank 3 is required as the cross-analysis comparison anchor")
    selected_indices, mad = select_top_mad(matrix, gene_ids, min(args.top_genes, len(gene_ids)))
    selected_gene_ids = tuple(gene_ids[index] for index in selected_indices)
    selected = matrix[selected_indices, :]

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent))
    try:
        write_tsv(
            stage / "selected_genes.tsv",
            [
                {
                    "selection_rank": rank,
                    "gene_id": gene_ids[index],
                    "gene_name": gene_names[index],
                    "mad": format(float(mad[index]), ".12g"),
                }
                for rank, index in enumerate(selected_indices, start=1)
            ],
        )
        results: dict[int, ConsensusNmfResult] = {}
        summary: list[dict[str, object]] = []
        studies = tuple(row["study_accession"] for row in metadata)
        for rank in ranks:
            result = run_consensus_nmf(
                selected,
                samples,
                rank=rank,
                runs=args.starts,
                seed_base=args.seed_base,
                max_iter=args.max_iter,
                tolerance=args.tolerance,
            )
            results[rank] = result
            _rank_outputs(stage, result, samples, metadata, selected_gene_ids)
            best = min(result.runs, key=lambda run: run.reconstruction_error)
            cluster_sizes = [int(np.sum(result.labels == state)) for state in range(rank)]
            association = _study_association(result.labels, studies, rank)
            summary.append(
                {
                    "rank": rank,
                    "sample_count": len(samples),
                    "gene_count": len(selected_gene_ids),
                    "nmf_runs": len(result.runs),
                    "best_fit": format(best.fit, ".12g"),
                    "median_fit": format(float(np.median([run.fit for run in result.runs])), ".12g"),
                    "cophenetic_correlation": format(result.cophenetic_correlation, ".12g"),
                    "dispersion": format(result.dispersion, ".12g"),
                    "consensus_sharpness": format(result.consensus_sharpness, ".12g"),
                    "average_silhouette": format(result.average_silhouette, ".12g"),
                    "minimum_cluster_size": min(cluster_sizes),
                    "cluster_sizes": ";".join(map(str, cluster_sizes)),
                    "converged_runs": sum(run.converged for run in result.runs),
                    "study_cramers_v": format(float(association["cramers_v"]), ".12g"),
                    "study_association_p": format(float(association["p_value"]), ".12g"),
                }
            )
        write_tsv(stage / "rank_summary.tsv", summary)
        eligible = [row for row in summary if int(row["minimum_cluster_size"]) >= 2]
        recommended = max(
            eligible or summary,
            key=lambda row: (
                float(row["cophenetic_correlation"]),
                float(row["average_silhouette"]),
                -int(row["rank"]),
            ),
        )
        counts_by_study = Counter(studies)
        receipt = {
            "status": "completed",
            "completed_at_utc": datetime.now(UTC).isoformat(),
            "analysis_id": args.analysis_id,
            "repo_commit": _git_commit(),
            "sample_policy": "primary_cohort_eligible=true; HGSOC tumour only",
            "input_value": "log1p_tpm from one shared Salmon/Gencode reference pipeline",
            "inputs": inputs,
            "sample_count": len(samples),
            "patient_count": len({row["patient_id"] for row in metadata}),
            "samples_by_study": dict(sorted(counts_by_study.items())),
            "gene_count_before_selection": len(gene_ids),
            "gene_count_after_selection": len(selected_gene_ids),
            "variable_gene_method": "median absolute deviation within this analysis matrix",
            "ranks": list(ranks),
            "random_initializations_per_rank": args.starts,
            "fixed_comparison_rank": 3,
            "recommended_rank": int(recommended["rank"]),
            "recommended_rank_rule": "highest cophenetic correlation among ranks with minimum cluster size >=2; tie by silhouette then lower rank",
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "scikit_learn": sklearn.__version__,
            },
            "claim_limit": "Response-blind technical NMF states; not Barnes subtype identities. Pooled results require explicit study-association audit.",
        }
        write_json(stage / "nmf_receipt.json", receipt)
        os.replace(stage, args.output_dir)
    except Exception as error:
        raise RuntimeError(f"NMF staging retained for inspection: {stage}") from error
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
