#!/usr/bin/env python3
"""Run a fail-closed Human-GEM WT versus ΔTPI1 FBA/FVA benchmark.

This consumes the preflight receipt rather than guessing a gene-to-reaction
mapping.  The model medium and objective are explicitly recorded.  It is a
model benchmark, not a replacement for Meeson's experimental validation, and
does not use any drug-response information.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Any

import cobra
from cobra.flux_analysis import flux_variability_analysis


def _finite(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} is non-numeric") from error
    if not math.isfinite(number):
        raise ValueError(f"{label} is non-finite")
    return number


def _preflight(path: Path) -> dict[str, Any]:
    root = json.loads(path.read_text(encoding="utf-8"))
    if root.get("status") != "valid" or root.get("response_blind") is not True:
        raise ValueError("preflight is not a valid response-blind receipt")
    deletion = root.get("gene_deletion")
    planned = root.get("planned_fva")
    model = root.get("model")
    if not isinstance(deletion, dict) or not isinstance(planned, dict) or not isinstance(model, dict):
        raise ValueError("preflight is structurally invalid")
    disabled = deletion.get("disabled_reaction_ids")
    targets = planned.get("reaction_ids")
    if not isinstance(disabled, list) or not disabled or not isinstance(targets, list) or not targets:
        raise ValueError("preflight has no deletion reactions or FVA targets")
    return root


def _solve(model: Any, label: str) -> dict[str, Any]:
    solution = model.optimize()
    status = str(solution.status)
    objective = _finite(solution.objective_value or 0.0, f"{label} objective")
    return {"status": status, "objective_value": objective}


def _ranges(model: Any, reactions: list[str], fraction: float) -> dict[str, list[float]]:
    if not 0 <= fraction <= 1:
        raise ValueError("fva fraction must be in [0, 1]")
    frame = flux_variability_analysis(
        model, reaction_list=reactions, fraction_of_optimum=fraction, processes=1
    )
    return {
        str(reaction): [_finite(row["minimum"], f"{reaction} minimum"), _finite(row["maximum"], f"{reaction} maximum")]
        for reaction, row in frame.iterrows()
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    preflight = _preflight(args.preflight)
    deleted = sorted(set(preflight["gene_deletion"]["disabled_reaction_ids"]))
    targets = sorted(set(preflight["planned_fva"]["reaction_ids"]))
    model = cobra.io.read_sbml_model(str(args.sbml))
    all_reactions = {reaction.id for reaction in model.reactions}
    absent = sorted((set(deleted) | set(targets) | {args.biomass_reaction}) - all_reactions)
    if absent:
        raise ValueError(f"model is missing planned reactions: {absent[:10]}")
    model.objective = args.biomass_reaction
    wild_type = _solve(model, "wild_type")
    if wild_type["status"].casefold() != "optimal" or wild_type["objective_value"] <= 0:
        raise ValueError("wild-type model has non-positive/non-optimal biomass; medium is not usable")
    wt_fva = _ranges(model, targets, args.fraction_of_optimum)

    deletion = model.copy()
    for reaction_id in deleted:
        deletion.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)
    delta = _solve(deletion, "delta_tpi1")
    if delta["status"].casefold() != "optimal":
        raise ValueError("delta-TPI1 model is not optimal")
    delta_fraction = args.fraction_of_optimum if delta["objective_value"] > 0 else 0.0
    delta_fva = _ranges(deletion, targets, delta_fraction)
    return {
        "status": "completed",
        "response_blind": True,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "preflight": str(args.preflight),
        "model": {"path": str(args.sbml), "biomass_reaction": args.biomass_reaction, "medium": model.medium},
        "gene_deletion": {
            "symbol": preflight["gene_deletion"]["symbol"],
            "disabled_reaction_ids": deleted,
            "implementation": "set GPR-validated disabled reaction bounds to zero",
        },
        "fba": {"wild_type": wild_type, "delta_tpi1": delta},
        "fva": {
            "target_reaction_count": len(targets),
            "wild_type_fraction_of_optimum": args.fraction_of_optimum,
            "delta_tpi1_fraction_of_optimum": delta_fraction,
            "wild_type_ranges": wt_fva,
            "delta_tpi1_ranges": delta_fva,
        },
        "claim_limit": (
            "Human-GEM model benchmark using the model's recorded medium; it is not an OCM-specific "
            "gene-dependency prediction, not a quantitative comparison to the published TPI1 table, "
            "and not a drug-response analysis."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--sbml", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--biomass-reaction", default="biomass_human")
    parser.add_argument("--fraction-of-optimum", type=float, default=0.9)
    args = parser.parse_args()
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"status": result["status"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
