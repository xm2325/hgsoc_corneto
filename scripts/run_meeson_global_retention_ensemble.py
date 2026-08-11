#!/usr/bin/env python3
"""Enumerate global-retention alternative optima for a reconstructed OCM cohort."""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

import cobra

from hgsoc_corneto.metabolic.global_retention import solve_global_retention
from hgsoc_corneto.metabolic.sequential import generate_meeson_candidates
from run_meeson_sequential_order_sensitivity import _expression, _read_receipt, _run_ids


def run(args: argparse.Namespace) -> dict[str, Any]:
    receipt = _read_receipt(args.receipt, args.study)
    run_ids = _run_ids(args.manifest, args.study)
    if set(run_ids) != set(receipt["targets"]):
        raise ValueError("receipt growth targets do not match frozen primary run IDs")
    if args.max_samples is not None:
        if args.max_samples < 1:
            raise ValueError("max-samples must be positive")
        run_ids = run_ids[: args.max_samples]
    expression = _expression(args.expression, run_ids, args.expression_transform)
    model = cobra.io.read_sbml_model(str(args.human_gem))
    results: list[dict[str, Any]] = []
    for run_id in run_ids:
        generated, _ = generate_meeson_candidates(
            model, expression[run_id].to_dict(), media_reactions=model.medium
        )
        candidates = [candidate for candidate in generated if candidate.reaction_id in receipt["selected"]]
        if {candidate.reaction_id for candidate in candidates} != receipt["selected"]:
            raise ValueError(f"{run_id}: receipt candidates cannot be regenerated")
        ensemble = solve_global_retention(
            model.copy(), candidates, biomass_id=args.biomass_reaction,
            growth_threshold=receipt["targets"][run_id], strict_margin=args.strict_margin,
            solver=args.solver, enumerate_alternatives=True, max_alternatives=args.max_alternatives,
        )
        alternatives = [list(solution.retained_reactions) for solution in ensemble.alternative_optima]
        frequency = Counter(reaction for solution in alternatives for reaction in solution)
        results.append(
            {
                "condition": run_id,
                "growth_threshold": receipt["targets"][run_id],
                "optimal_retained_count": ensemble.optimal_retained_count,
                "solution_count": len(alternatives),
                "alternatives_truncated": ensemble.alternatives_truncated,
                "retained_sets": alternatives,
                "retention_frequency": {
                    reaction: frequency[reaction] / len(alternatives)
                    for reaction in sorted(receipt["selected"])
                },
            }
        )
    return {
        "status": "completed",
        "response_blind": True,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "study_accession": args.study,
        "sample_count": len(results),
        "formulation": "CORNETO simultaneous constraint-retention MILP with no-good-cut optimum enumeration",
        "max_alternatives": args.max_alternatives,
        "strict_margin": args.strict_margin,
        "solver": args.solver,
        "receipt": str(args.receipt),
        "conditions": results,
        "claim_limit": (
            "Alternative-optimum ensemble for the new global-retention benchmark; it is not an exact "
            "Meeson MitoCore/sFBA reproduction and does not use response data."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--biomass-reaction", default="biomass_human")
    parser.add_argument("--expression-transform", choices=("raw_tpm", "log1p_tpm"), default="log1p_tpm")
    parser.add_argument("--max-alternatives", type=int, default=30)
    parser.add_argument("--max-samples", type=int)
    parser.add_argument("--strict-margin", type=float, default=1e-6)
    parser.add_argument("--solver", default="gurobi")
    args = parser.parse_args()
    if args.max_alternatives < 1:
        parser.error("max-alternatives must be positive")
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(args.output.suffix + ".tmp")
    tmp.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(args.output)
    print(json.dumps({"status": result["status"], "study": args.study, "output": str(args.output)}))


if __name__ == "__main__":
    main()
