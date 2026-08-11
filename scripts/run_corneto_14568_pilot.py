#!/usr/bin/env python3
"""Run a response-blind CORNETO/Human-GEM pilot on E-MTAB-14568.

This is a solver/data-integration pilot, not a paclitaxel-response analysis.
The repository currently has no real-data regulatory-network implementation;
this entry point therefore exercises the pinned metabolic CORNETO formulation
against the harmonised 14568 TPM matrix and records every non-default choice.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any

import cobra
import pandas as pd

from hgsoc_corneto.metabolic.joint_fba import compare_independent_and_joint_sparse_fba
from hgsoc_corneto.metabolic.sequential import generate_meeson_candidates


def _read_manifest(path: Path, *, study: str, primary_only: bool) -> list[dict[str, str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame = frame[frame["study_accession"] == study]
    if primary_only:
        frame = frame[
            (frame["sample_class"] == "tumour")
            & (frame["primary_cohort_eligible"] == "true")
        ]
    rows = frame.to_dict(orient="records")
    rows.sort(key=lambda row: row["run_accession"])
    if not rows:
        raise ValueError("No runs remain after manifest filtering")
    return rows


def _read_expression(path: Path, run_ids: list[str], transform: str) -> pd.DataFrame:
    usecols = ["gene_id", "gene_name", *run_ids]
    frame = pd.read_csv(path, sep="\t", compression="gzip", usecols=usecols)
    frame["gene_id"] = frame["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    frame = frame.drop(columns=["gene_name"]).set_index("gene_id")
    frame = frame.apply(pd.to_numeric, errors="raise")
    if frame.index.duplicated().any():
        frame = frame.groupby(level=0, sort=False).sum()
    if transform == "raw_tpm":
        pass
    elif transform == "log1p_tpm":
        frame = frame.map(lambda value: math.log1p(float(value)))
    elif transform == "log1p_tpm_div5":
        frame = frame.map(lambda value: math.log1p(float(value)) / 5.0)
    else:  # pragma: no cover - argparse protects this branch
        raise ValueError(f"Unknown expression transform: {transform}")
    if not frame.index.is_unique:
        raise ValueError("Gene IDs are not unique after version stripping")
    return frame


def _solver_choice(requested: str) -> tuple[str, list[str], str | None]:
    import corneto as cn

    available = [str(value) for value in cn.opt.available_solvers()]
    normalized = {value.upper() for value in available}
    if requested == "auto":
        if "GUROBI" in normalized and importlib.util.find_spec("gurobipy") is not None:
            return "gurobi", available, None
        if "MOSEK" in normalized and importlib.util.find_spec("mosek") is not None:
            return "mosek", available, None
        return "highs", available, "Commercial solver package/license not visible; using HiGHS"
    if requested.upper() not in normalized:
        raise RuntimeError(
            f"Requested solver {requested!r} is unavailable; available={available}"
        )
    return requested, available, None


def _candidate_sets(
    model: Any,
    expression: pd.DataFrame,
    run_ids: list[str],
    *,
    max_candidates: int,
) -> tuple[
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
    list[str],
    dict[str, float],
]:
    per_sample: dict[str, dict[str, Any]] = {}
    audits: dict[str, dict[str, Any]] = {}
    for run_id in run_ids:
        candidates, audit = generate_meeson_candidates(
            model,
            expression[run_id].to_dict(),
            media_reactions=model.medium,
        )
        positive = {
            candidate.reaction_id: candidate
            for candidate in candidates
            if candidate.proposed_upper > 0
        }
        per_sample[run_id] = positive
        audits[run_id] = {
            "all_candidates": len(candidates),
            "positive_candidates": len(positive),
            "category_counts": audit["candidate_counts_by_category"],
            "skipped_counts": audit["skipped_counts"],
        }

    scores: list[tuple[float, str]] = []
    for reaction_id in sorted({rid for values in per_sample.values() for rid in values}):
        bounds = [
            float(per_sample[run_id][reaction_id].proposed_upper)
            for run_id in run_ids
            if reaction_id in per_sample[run_id]
        ]
        scores.append((statistics.median(bounds), reaction_id))
    scores.sort()
    selected = [reaction_id for _, reaction_id in scores[:max_candidates]]
    if not selected:
        raise ValueError("No positive expression-derived reaction candidates remain")
    selected_scores = {reaction_id: score for score, reaction_id in scores[:max_candidates]}
    return per_sample, audits, selected, selected_scores


def _reaction_bounds(
    model: Any,
    per_sample: dict[str, dict[str, Any]],
    run_ids: list[str],
    selected: list[str],
    growth_fraction: float,
) -> tuple[
    dict[str, dict[str, tuple[float, float]]],
    dict[str, float],
    dict[str, float],
    dict[str, int],
]:
    bounds: dict[str, dict[str, tuple[float, float]]] = {}
    optima: dict[str, float] = {}
    targets: dict[str, float] = {}
    skipped: dict[str, int] = {}
    for run_id in run_ids:
        sample_bounds: dict[str, tuple[float, float]] = {}
        skipped_count = 0
        for reaction_id in selected:
            candidate = per_sample[run_id].get(reaction_id)
            if candidate is None:
                continue
            reaction = model.reactions.get_by_id(reaction_id)
            lower = max(float(reaction.lower_bound), float(candidate.proposed_lower))
            upper = min(float(reaction.upper_bound), float(candidate.proposed_upper))
            if lower > upper:
                skipped_count += 1
                continue
            sample_bounds[reaction_id] = (lower, upper)

        constrained = model.copy()
        for reaction_id, reaction_bounds in sample_bounds.items():
            constrained.reactions.get_by_id(reaction_id).bounds = reaction_bounds
        constrained.objective = "biomass_human"
        optimum_solution = constrained.optimize()
        optimum = float(optimum_solution.objective_value or 0.0)
        if not math.isfinite(optimum) or optimum <= 0:
            raise RuntimeError(
                f"Non-positive or non-finite biomass optimum for {run_id}: {optimum}"
            )
        biomass = model.reactions.get_by_id("biomass_human")
        target = growth_fraction * optimum
        sample_bounds["biomass_human"] = (
            max(float(biomass.lower_bound), target),
            float(biomass.upper_bound),
        )
        bounds[run_id] = sample_bounds
        optima[run_id] = optimum
        targets[run_id] = target
        skipped[run_id] = skipped_count
    return bounds, optima, targets, skipped


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--study", default="E-MTAB-14568")
    parser.add_argument("--primary-only", action="store_true")
    parser.add_argument("--max-samples", type=int, default=3)
    parser.add_argument("--max-candidates", type=int, default=250)
    parser.add_argument(
        "--growth-fraction",
        type=float,
        default=0.9,
        help="Minimum fraction of each expression-constrained biomass optimum.",
    )
    parser.add_argument(
        "--expression-transform",
        choices=("raw_tpm", "log1p_tpm", "log1p_tpm_div5"),
        default="log1p_tpm",
    )
    parser.add_argument("--solver", choices=("auto", "highs", "mosek", "gurobi"), default="auto")
    parser.add_argument("--independent-lambda", type=float, default=0.1)
    parser.add_argument("--joint-lambda", type=float, default=1.0)
    args = parser.parse_args()

    rows = _read_manifest(args.manifest, study=args.study, primary_only=args.primary_only)
    if args.max_samples < 1:
        raise ValueError("--max-samples must be positive")
    rows = rows[: args.max_samples]
    run_ids = [row["run_accession"] for row in rows]
    expression = _read_expression(args.expression, run_ids, args.expression_transform)
    model = cobra.io.read_sbml_model(str(args.human_gem))
    requested_solver = args.solver
    solver, available_solvers, fallback_reason = _solver_choice(requested_solver)

    per_sample, candidate_audits, selected, selected_scores = _candidate_sets(
        model, expression, run_ids, max_candidates=args.max_candidates
    )
    if not 0 < args.growth_fraction <= 1:
        raise ValueError("--growth-fraction must be in (0, 1]")
    reaction_bounds, growth_optima, growth_targets, bounds_skipped = _reaction_bounds(
        model, per_sample, run_ids, selected, args.growth_fraction
    )
    objectives = {run_id: {"biomass_human": -1.0} for run_id in run_ids}

    solver_fallback: str | None = None
    try:
        comparison = compare_independent_and_joint_sparse_fba(
            model,
            objectives=objectives,
            reaction_bounds=reaction_bounds,
            independent_lambda=args.independent_lambda,
            joint_lambda=args.joint_lambda,
            solver=solver,
        )
    except Exception as error:
        error_text = str(error).casefold()
        if (
            requested_solver == "auto"
            and solver == "gurobi"
            and ("gurobi" in error_text or "license" in error_text)
        ):
            solver_fallback = f"Gurobi solve failed: {type(error).__name__}"
            solver = "highs"
            comparison = compare_independent_and_joint_sparse_fba(
                model,
                objectives=objectives,
                reaction_bounds=reaction_bounds,
                independent_lambda=args.independent_lambda,
                joint_lambda=args.joint_lambda,
                solver=solver,
            )
        else:
            raise

    result = {
        "status": "completed",
        "study_accession": args.study,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "repo_commit": os.environ.get("HGSOC_CORNETO_REPO_COMMIT"),
        "primary_only": args.primary_only,
        "sample_count": len(run_ids),
        "samples": rows,
        "expression": {
            "path": str(args.expression),
            "transform": args.expression_transform,
            "gene_count": int(expression.shape[0]),
        },
        "model": {
            "path": str(args.human_gem),
            "reactions": len(model.reactions),
            "genes": len(model.genes),
            "biomass_id": "biomass_human",
        },
        "solver": {
            "requested": requested_solver,
            "used": solver,
            "available": available_solvers,
            "fallback_reason": fallback_reason,
            "solve_fallback": solver_fallback,
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
            "independent_lambda": args.independent_lambda,
            "joint_lambda": args.joint_lambda,
        },
        "corneto": comparison.to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
