#!/usr/bin/env python3
"""Run joint-only metabolic CORNETO on the frozen pooled primary cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import cobra
import pandas as pd

from hgsoc_corneto.metabolic.joint_fba import solve_joint_sparse_fba
from run_corneto_14568_pilot import (
    _candidate_sets,
    _reaction_bounds,
    _read_expression,
    _solver_choice,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_samples(path: Path, policy: str, expected_samples: int) -> list[dict[str, str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "study_accession",
        "run_accession",
        "canonical_ocm_id",
        "patient_id",
        "sample_class",
        "histotype_group",
        "primary_cohort_eligible",
    }
    if missing := required - set(frame.columns):
        raise ValueError(f"pooled sample manifest is missing fields: {sorted(missing)}")
    if policy == "one_per_study":
        frame = frame.groupby("study_accession", sort=False, as_index=False).head(1)
    elif policy != "all":
        raise ValueError(f"unknown sample policy {policy!r}")
    rows: list[dict[str, str]] = frame.to_dict(orient="records")
    if len(rows) != expected_samples:
        raise ValueError(f"expected {expected_samples} selected samples, found {len(rows)}")
    if any(
        row["sample_class"] != "tumour"
        or row["histotype_group"] != "HGSOC"
        or row["primary_cohort_eligible"].casefold() != "true"
        for row in rows
    ):
        raise ValueError("pooled selection contains a non-primary HGSOC tumour row")
    run_ids = [row["run_accession"] for row in rows]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("pooled selection contains duplicate run accessions")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--pooled-receipt", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-policy", choices=("all", "one_per_study"), required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--max-candidates", type=int, default=25)
    parser.add_argument("--growth-fraction", type=float, default=0.9)
    parser.add_argument("--joint-lambda", type=float, default=1.0)
    parser.add_argument("--solver", choices=("gurobi", "highs", "mosek"), default="gurobi")
    args = parser.parse_args()

    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite output: {args.output}")
    if args.max_candidates < 1:
        raise ValueError("--max-candidates must be positive")
    if not 0 < args.growth_fraction <= 1:
        raise ValueError("--growth-fraction must be in (0, 1]")
    if args.joint_lambda < 0:
        raise ValueError("--joint-lambda must be non-negative")
    pooled_receipt: Any = json.loads(args.pooled_receipt.read_text(encoding="utf-8"))
    if not isinstance(pooled_receipt, dict) or pooled_receipt.get("status") != "completed":
        raise ValueError("pooled expression receipt is not completed")
    rows = _read_samples(args.sample_manifest, args.sample_policy, args.expected_samples)
    run_ids = [row["run_accession"] for row in rows]
    expression = _read_expression(args.expression, run_ids, "log1p_tpm")
    expression = expression.loc[:, run_ids]
    model = cobra.io.read_sbml_model(str(args.human_gem))
    solver, available_solvers, fallback_reason = _solver_choice(args.solver)
    per_sample, candidate_audits, selected, selected_scores = _candidate_sets(
        model,
        expression,
        run_ids,
        max_candidates=args.max_candidates,
    )
    reaction_bounds, growth_optima, growth_targets, bounds_skipped = _reaction_bounds(
        model,
        per_sample,
        run_ids,
        selected,
        args.growth_fraction,
    )
    objectives = {run_id: {"biomass_human": -1.0} for run_id in run_ids}
    joint = solve_joint_sparse_fba(
        model,
        objectives=objectives,
        reaction_bounds=reaction_bounds,
        joint_lambda=args.joint_lambda,
        solver=solver,
    )
    result = {
        "status": "completed",
        "analysis": "pooled_primary_joint_metabolic_sparse_fba",
        "joint_only": True,
        "response_blind": True,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "repo_commit": os.environ.get("HGSOC_CORNETO_REPO_COMMIT"),
        "sample_policy": args.sample_policy,
        "sample_count": len(run_ids),
        "samples": rows,
        "expression": {
            "path": str(args.expression),
            "transform": "log1p_tpm",
            "gene_count": int(expression.shape[0]),
        },
        "provenance": {
            "expression_sha256": _sha256(args.expression),
            "sample_manifest_sha256": _sha256(args.sample_manifest),
            "pooled_receipt_sha256": _sha256(args.pooled_receipt),
            "human_gem_sha256": _sha256(args.human_gem),
            "runner_sha256": _sha256(Path(__file__)),
        },
        "model": {
            "path": str(args.human_gem),
            "reactions": len(model.reactions),
            "genes": len(model.genes),
            "biomass_id": "biomass_human",
        },
        "solver": {
            "requested": args.solver,
            "used": solver,
            "available": available_solvers,
            "fallback_reason": fallback_reason,
        },
        "candidate_selection": {
            "policy": "positive candidates ranked by median proposed upper bound; lowest retained",
            "max_candidates": args.max_candidates,
            "selected_count": len(selected),
            "selected_reaction_ids": selected,
            "selected_median_proposed_upper": selected_scores,
            "audits_by_sample": candidate_audits,
            "bounds_clamped_to_human_gem": True,
            "bounds_skipped_as_empty_intersection": bounds_skipped,
        },
        "objective": {
            "biomass_coefficient": -1.0,
            "growth_fraction": args.growth_fraction,
            "growth_optima": growth_optima,
            "growth_targets": growth_targets,
            "joint_lambda": args.joint_lambda,
        },
        "corneto": joint.to_dict(),
        "claim_limit": (
            "Pooled response-blind model-predicted feasible metabolic state; not measured "
            "flux, phenotype association, or evidence that pooling is superior to cohort models."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"status": "completed", "samples": len(run_ids), "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
