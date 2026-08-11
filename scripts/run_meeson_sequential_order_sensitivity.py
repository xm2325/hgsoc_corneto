#!/usr/bin/env python3
"""Test order sensitivity of the public Meeson sequential rule on an OCM cohort.

The analysis reconstructs candidate bounds from the frozen expression matrix,
uses the candidate set and growth targets recorded by a completed metabolic
receipt, and applies the published sequential rule under a canonical and
seeded random reaction order.  It is response-blind.  It is not an exact
reproduction of Meeson's 49-OCM models because their expression-plus-growth
input is not public.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
from pathlib import Path
from typing import Any

import cobra
import pandas as pd

from hgsoc_corneto.metabolic.sequential import (
    CandidateConstraint,
    apply_sequential_constraints,
    generate_meeson_candidates,
)


def _read_receipt(path: Path, study: str) -> dict[str, Any]:
    root = json.loads(path.read_text(encoding="utf-8"))
    if root.get("status") != "completed" or root.get("study_accession") != study:
        raise ValueError("receipt is not a completed matching cohort result")
    selected = root.get("candidate_selection", {}).get("selected_reaction_ids")
    targets = root.get("objective", {}).get("growth_targets")
    if not isinstance(selected, list) or not selected or len(selected) != len(set(selected)):
        raise ValueError("receipt has invalid selected_reaction_ids")
    if not isinstance(targets, dict):
        raise ValueError("receipt has invalid growth_targets")
    return {"selected": set(selected), "targets": {str(key): float(value) for key, value in targets.items()}}


def _run_ids(manifest: Path, study: str) -> list[str]:
    frame = pd.read_csv(manifest, sep="\t", dtype=str, keep_default_na=False)
    subset = frame[
        (frame["study_accession"] == study)
        & (frame["sample_class"] == "tumour")
        & (frame["histotype_group"] == "HGSOC")
        & (frame["primary_cohort_eligible"] == "true")
        & (frame["is_representative_rna_library"] == "true")
    ]
    values = sorted(subset["run_accession"].tolist())
    if not values or len(values) != len(set(values)):
        raise ValueError("primary run IDs are missing or duplicated")
    return values


def _expression(path: Path, run_ids: list[str], transform: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="gzip", usecols=["gene_id", "gene_name", *run_ids])
    frame["gene_id"] = frame["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    frame = frame.drop(columns=["gene_name"]).set_index("gene_id").apply(pd.to_numeric, errors="raise")
    if frame.index.duplicated().any():
        frame = frame.groupby(level=0, sort=False).sum()
    if transform == "log1p_tpm":
        frame = frame.map(lambda value: math.log1p(float(value)))
    elif transform != "raw_tpm":
        raise ValueError(f"unsupported transform {transform!r}")
    return frame


def _jaccard(left: set[str], right: set[str]) -> float | None:
    union = left | right
    return None if not union else len(left & right) / len(union)


def _summary(result: Any) -> dict[str, Any]:
    retained = set(result.retained_reactions)
    return {
        "retained_reaction_ids": sorted(retained),
        "retained_count": len(retained),
        "final_status": result.final_status,
        "final_growth": result.final_growth,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    receipt = _read_receipt(args.receipt, args.study)
    run_ids = _run_ids(args.manifest, args.study)
    if set(run_ids) != set(receipt["targets"]):
        raise ValueError("receipt growth targets do not match frozen primary run IDs")
    expression = _expression(args.expression, run_ids, args.expression_transform)
    model = cobra.io.read_sbml_model(str(args.human_gem))
    if args.biomass_reaction not in {reaction.id for reaction in model.reactions}:
        raise ValueError("biomass reaction is absent from model")
    if args.permutations < 1:
        raise ValueError("permutations must be positive")

    conditions: list[dict[str, Any]] = []
    for condition_index, run_id in enumerate(run_ids):
        generated, _ = generate_meeson_candidates(
            model, expression[run_id].to_dict(), media_reactions=model.medium
        )
        candidates = [candidate for candidate in generated if candidate.reaction_id in receipt["selected"]]
        candidate_ids = {candidate.reaction_id for candidate in candidates}
        if candidate_ids != receipt["selected"]:
            absent = sorted(receipt["selected"] - candidate_ids)
            raise ValueError(f"{run_id}: receipt candidates cannot be regenerated: {absent[:8]}")
        canonical = apply_sequential_constraints(
            model.copy(), candidates, biomass_id=args.biomass_reaction,
            growth_threshold=receipt["targets"][run_id], semantics=args.semantics,
        )
        canonical_summary = _summary(canonical)
        alternatives: list[dict[str, Any]] = []
        for permutation in range(args.permutations):
            shuffled = list(candidates)
            random.Random(args.seed + condition_index * 100_003 + permutation).shuffle(shuffled)
            alternative = apply_sequential_constraints(
                model.copy(), shuffled, biomass_id=args.biomass_reaction,
                growth_threshold=receipt["targets"][run_id], semantics=args.semantics,
            )
            item = _summary(alternative)
            item["permutation"] = permutation
            item["retained_jaccard_to_canonical"] = _jaccard(
                set(canonical_summary["retained_reaction_ids"]), set(item["retained_reaction_ids"])
            )
            alternatives.append(item)
        jaccards = [item["retained_jaccard_to_canonical"] for item in alternatives]
        counts = [item["retained_count"] for item in alternatives]
        conditions.append(
            {
                "condition": run_id,
                "growth_threshold": receipt["targets"][run_id],
                "candidate_count": len(candidates),
                "canonical": canonical_summary,
                "permutations": alternatives,
                "summary": {
                    "retained_count_min": min(counts),
                    "retained_count_max": max(counts),
                    "retained_count_unique": sorted(set(counts)),
                    "jaccard_to_canonical_min": min(jaccards),
                    "jaccard_to_canonical_mean": sum(jaccards) / len(jaccards),
                    "jaccard_to_canonical_max": max(jaccards),
                },
            }
        )
    return {
        "status": "completed",
        "response_blind": True,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "study_accession": args.study,
        "sample_count": len(run_ids),
        "semantics": args.semantics,
        "permutations_per_sample": args.permutations,
        "random_seed": args.seed,
        "receipt": str(args.receipt),
        "claim_limit": (
            "Order-sensitivity benchmark of a public sequential rule on reconstructed OCM inputs; "
            "not an exact reproduction of Meeson's unavailable 49-OCM inputs and not a drug-response analysis."
        ),
        "conditions": conditions,
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
    parser.add_argument("--semantics", choices=("published", "bounds_safe"), default="bounds_safe")
    parser.add_argument("--permutations", type=int, default=25)
    parser.add_argument("--seed", type=int, default=20260811)
    args = parser.parse_args()
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(args.output.suffix + ".tmp")
    tmp.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(args.output)
    print(json.dumps({"status": result["status"], "study": args.study, "output": str(args.output)}))


if __name__ == "__main__":
    main()
