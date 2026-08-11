#!/usr/bin/env python3
"""Compare pooled-60 rank-3 NMF with primary-only cohort-stratified rank-3 fits."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import chi2_contingency, spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


STUDIES = ("E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_w(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    return [row[0] for row in rows], header[1:], np.asarray([[float(x) for x in row[1:]] for row in rows])


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError(f"refusing to overwrite: {args.output_dir}")
    pooled_dir = args.repo / "data/processed/rna/pooled_primary60/nmf/rank_3"
    pooled_assign = {row["run_accession"]: row for row in read_tsv(pooled_dir / "assignments.tsv")}
    p_genes, p_states, p_w = load_w(pooled_dir / "best_w_gene_loadings.tsv.gz")
    p_index = {gene: i for i, gene in enumerate(p_genes)}
    mappings: list[dict[str, object]] = []
    metrics: list[dict[str, object]] = []
    for study in STUDIES:
        cohort_dir = args.repo / "data/processed/rna" / study / "nmf_primary/rank_3"
        cohort_assign = {row["run_accession"]: row for row in read_tsv(cohort_dir / "assignments.tsv")}
        c_genes, c_states, c_w = load_w(cohort_dir / "best_w_gene_loadings.tsv.gz")
        common = sorted(set(p_genes).intersection(c_genes))
        if len(common) < 1000:
            raise ValueError(f"too few common loading genes for {study}: {len(common)}")
        pi = [p_index[gene] for gene in common]
        ci_map = {gene: i for i, gene in enumerate(c_genes)}
        ci = [ci_map[gene] for gene in common]
        correlation = np.empty((len(c_states), len(p_states)))
        for i in range(len(c_states)):
            for j in range(len(p_states)):
                correlation[i, j] = float(spearmanr(c_w[ci, i], p_w[pi, j]).statistic)
        row_ind, col_ind = linear_sum_assignment(-correlation)
        state_map = {c_states[i]: p_states[j] for i, j in zip(row_ind, col_ind, strict=True)}
        for i, j in zip(row_ind, col_ind, strict=True):
            mappings.append({
                "study_accession": study,
                "cohort_state": c_states[i],
                "pooled_state": p_states[j],
                "loading_spearman": format(float(correlation[i, j]), ".12g"),
                "common_genes": len(common),
            })
        runs = sorted(cohort_assign)
        if any(run not in pooled_assign for run in runs):
            raise ValueError(f"cohort run absent from pooled assignments: {study}")
        cohort_labels = [cohort_assign[run]["state"] for run in runs]
        pooled_labels = [pooled_assign[run]["state"] for run in runs]
        mapped_labels = [state_map[label] for label in cohort_labels]
        metrics.append({
            "study_accession": study,
            "sample_count": len(runs),
            "adjusted_rand_index": format(float(adjusted_rand_score(cohort_labels, pooled_labels)), ".12g"),
            "normalized_mutual_information": format(float(normalized_mutual_info_score(cohort_labels, pooled_labels)), ".12g"),
            "mapped_assignment_agreement": format(float(np.mean(np.asarray(mapped_labels) == np.asarray(pooled_labels))), ".12g"),
            "mean_matched_loading_spearman": format(float(np.mean([correlation[i, j] for i, j in zip(row_ind, col_ind, strict=True)])), ".12g"),
            "minimum_matched_loading_spearman": format(float(np.min([correlation[i, j] for i, j in zip(row_ind, col_ind, strict=True)])), ".12g"),
            "common_loading_genes": len(common),
        })
    study_levels = list(STUDIES)
    state_levels = sorted({row["state"] for row in pooled_assign.values()})
    contingency = np.zeros((len(study_levels), len(state_levels)), dtype=int)
    for row in pooled_assign.values():
        contingency[study_levels.index(row["study_accession"]), state_levels.index(row["state"])] += 1
    chi2, pvalue, _, _ = chi2_contingency(contingency)
    n = int(contingency.sum())
    cramers_v = float(np.sqrt(chi2 / (n * min(contingency.shape[0] - 1, contingency.shape[1] - 1))))
    args.output_dir.mkdir(parents=True)
    write_tsv(args.output_dir / "state_mapping.tsv", mappings)
    write_tsv(args.output_dir / "cohort_comparison.tsv", metrics)
    receipt = {
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "rank": 3,
        "pooled_sample_count": len(pooled_assign),
        "cohort_metrics": metrics,
        "pooled_state_by_study": {
            "study_levels": study_levels,
            "state_levels": state_levels,
            "counts": contingency.tolist(),
            "chi_square": float(chi2),
            "p_value": float(pvalue),
            "cramers_v": cramers_v,
        },
        "interpretation_guardrail": "strong state-study association indicates pooled NMF may be dominated by accession/batch structure",
    }
    (args.output_dir / "comparison_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
