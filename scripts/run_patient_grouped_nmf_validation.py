#!/usr/bin/env python3
"""Train-only NMF feature selection/projection and repeated patient balancing.

No response labels, Gurobi, metabolic receipts or pooled-fit bases enter the
held-out fits. Rank 2 and 3 are prespecified sensitivity analyses, not a
test-selected winning rank. Projection quality alone is not subtype replication.
"""

from __future__ import annotations

import argparse
import platform
import warnings
from pathlib import Path

import numpy as np
import sklearn
from research_validation_common import (
    STUDIES,
    heldout_split,
    load_data,
    one_per_patient,
    primary_data,
    select_mad,
    sha,
    stamp,
    tsv,
    write_json,
)
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import adjusted_rand_score


def normalized_labels(w, h):
    # NMF factors are scale-indeterminate: classify contributions, not raw W.
    contribution = w * h.sum(axis=1)[None, :]
    return np.argmax(contribution, axis=1).astype(int)


def fit_best(x, rank, starts, seed, max_iter):
    best, diagnostics = None, []
    for start in range(starts):
        model = NMF(
            n_components=rank,
            init="random",
            random_state=seed + start,
            max_iter=max_iter,
            tol=1e-4,
            solver="cd",
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ConvergenceWarning)
            w = model.fit_transform(x)
        converged = not any(issubclass(v.category, ConvergenceWarning) for v in caught)
        diagnostics.append(
            {
                "seed": seed + start,
                "converged": converged,
                "iterations": int(model.n_iter_),
                "error": float(model.reconstruction_err_),
            }
        )
        if converged and (best is None or model.reconstruction_err_ < best[0].reconstruction_err_):
            best = model, w
    if best is None:
        raise RuntimeError("no converged training start; do not interpret projection")
    return best[0], best[1], diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--rna-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-genes", type=int, default=6000)
    parser.add_argument("--starts", type=int, default=20)
    parser.add_argument("--draws", type=int, default=30)
    parser.add_argument("--max-iter", type=int, default=3000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "patient_grouped_nmf_validation.v1",
        "status": "running",
        "started_at_utc": stamp(),
        "parameters": vars(args).copy(),
        "folds": [],
        "patient_balanced_draws": [],
        "gurobi_sessions_requested": 0,
        "script_sha256": sha(__file__),
        "common_code_sha256": sha(Path(__file__).with_name("research_validation_common.py")),
        "software": {"python": platform.python_version(), "sklearn": sklearn.__version__},
        "claim_limit": (
            "Response-blind transferability/robustness diagnostics, not clinical subtypes "
            "or independent biological validation."
        ),
    }
    report["parameters"] = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    output = args.output_dir / "receipt.json"
    write_json(output, report)
    try:
        genes, names, x, metadata, sources = load_data(args.root, args.rna_root)
        x, metadata = primary_data(x, metadata)
        report["sources"] = sources
        report["sample_count"], report["patient_count"] = (
            len(metadata),
            len({m["patient_id"] for m in metadata}),
        )
        for fold, study in enumerate(STUDIES):
            train, test = heldout_split(metadata, study, 20260908 + fold)
            selected = select_mad(x[:, train], genes, args.top_genes)
            a, b = x[selected][:, train].T, x[selected][:, test].T
            for rank in (2, 3):
                model, w, starts = fit_best(
                    a, rank, args.starts, 20260908 + 1000 * fold + rank * 100, args.max_iter
                )
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always", ConvergenceWarning)
                    projected = model.transform(b)
                if any(issubclass(v.category, ConvergenceWarning) for v in caught):
                    raise RuntimeError("held-out transform failed to converge")
                predicted = projected @ model.components_
                error = np.square(b - predicted).sum(axis=1)
                baseline = np.square(b - a.mean(axis=0)).sum(axis=1)
                train_patients = sorted({metadata[i]["patient_id"] for i in train})
                test_patients = sorted({metadata[i]["patient_id"] for i in test})
                per_sample = [
                    {
                        "run_accession": metadata[i]["run_accession"],
                        "patient_id": metadata[i]["patient_id"],
                        "projected_state": int(normalized_labels(projected, model.components_)[j])
                        + 1,
                        "squared_error": float(error[j]),
                        "training_mean_squared_error": float(baseline[j]),
                        "relative_improvement_over_training_mean": float(1 - error[j] / baseline[j])
                        if baseline[j] > 0
                        else None,
                    }
                    for j, i in enumerate(test)
                ]
                family_scores = [
                    np.mean(
                        [
                            r["relative_improvement_over_training_mean"]
                            for r in per_sample
                            if r["patient_id"] == p
                            and r["relative_improvement_over_training_mean"] is not None
                        ]
                    )
                    for p in test_patients
                ]
                if not np.isfinite(family_scores).all():
                    raise ValueError("undefined patient-level projection metric")
                fold_result = {
                    "heldout_study": study,
                    "rank": rank,
                    "training_patients": train_patients,
                    "heldout_patients": test_patients,
                    "patient_overlap": [],
                    "training_runs": [metadata[i]["run_accession"] for i in train],
                    "training_gene_ids": [genes[i] for i in selected],
                    "starts": starts,
                    "patient_mean_relative_improvement": float(np.mean(family_scores)),
                    "heldout_samples": per_sample,
                    "claim_limit": (
                        "Reconstruction vs training-mean baseline is necessary diagnostic "
                        "evidence, not proof of subtype replication."
                    ),
                }
                report["folds"].append(fold_result)
                np.savez_compressed(
                    args.output_dir / (study + "_rank" + str(rank) + ".npz"),
                    genes=np.asarray([genes[i] for i in selected]),
                    basis=model.components_,
                    heldout_contributions=projected * model.components_.sum(axis=1)[None, :],
                )
                write_json(output, report)
                print(study, rank, "projection recorded", flush=True)
        all_indices = list(range(len(metadata)))
        for rank in (2, 3):
            reference_path = (
                args.root
                / "data/processed/rna/pooled_primary60/nmf"
                / f"rank_{rank}"
                / "assignments.tsv"
            )
            reference = {r["run_accession"]: r["state"] for r in tsv(reference_path)}
            if set(reference) != {m["run_accession"] for m in metadata}:
                raise ValueError("historical assignment coverage mismatch")
            report.setdefault("reference_assignments", []).append(
                {"path": str(reference_path), "sha256": sha(reference_path)}
            )
            for draw in range(args.draws):
                keep = one_per_patient(metadata, all_indices, 20260908 + draw)
                selected = select_mad(x[:, keep], genes, args.top_genes)
                model, w, starts = fit_best(
                    x[selected][:, keep].T,
                    rank,
                    args.starts,
                    20260908 + 10000 * rank + draw * 100,
                    args.max_iter,
                )
                labels = normalized_labels(w, model.components_)
                ari = adjusted_rand_score(
                    [reference[metadata[i]["run_accession"]] for i in keep], labels
                )
                report["patient_balanced_draws"].append(
                    {
                        "rank": rank,
                        "draw": draw,
                        "seed": 20260908 + draw,
                        "runs": [metadata[i]["run_accession"] for i in keep],
                        "patients": [metadata[i]["patient_id"] for i in keep],
                        "ari_vs_historical_pooled": float(ari),
                        "starts": starts,
                        "claim_limit": (
                            "Same-data sensitivity; also changes feature selection and uses "
                            "scale-normalized factor contributions, not legacy consensus labels."
                        ),
                    }
                )
                write_json(output, report)
        report["status"] = "completed"
    except Exception as exc:
        report["status"], report["error"] = "incomplete", repr(exc)
        raise
    finally:
        report["finished_at_utc"] = stamp()
        write_json(output, report)
    print("completed", output, flush=True)


if __name__ == "__main__":
    main()
