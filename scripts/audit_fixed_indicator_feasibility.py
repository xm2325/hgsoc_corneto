#!/usr/bin/env python3
"""Test saved reaction selections with continuous HiGHS LPs (no Gurobi).

No original file is changed. A feasible restricted LP is not validation of
medium composition, expression interpretation or disease specificity.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from audit_metabolic_scientific_validity import STUDIES, read_json, read_sbml, sha256
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def restricted_lp(reactions, sample_bounds, selected, required_biomass=None):
    unknown = set(selected) - set(reactions)
    if unknown:
        raise ValueError(f"unknown selected reactions: {sorted(unknown)}")
    ids = sorted(set(selected))
    bounds = []
    effective = {}
    for rid, reaction in reactions.items():
        lower, upper = reaction["lower"], reaction["upper"]
        if rid in sample_bounds:
            lower = max(lower, sample_bounds[rid][0])
            upper = min(upper, sample_bounds[rid][1])
        if rid == "biomass_human" and required_biomass is not None:
            lower = max(lower, required_biomass)
        effective[rid] = (lower, upper)
        if lower > upper or (rid not in selected and (lower > 0 or upper < 0)):
            return {
                "status": "infeasible",
                "reason": "off_reaction_requires_nonzero_flux",
                "reaction": rid,
                "success": False,
            }
    metabolite_ids = sorted({m for r in ids for m in reactions[r]["stoichiometry"]})
    metabolite_index = {m: i for i, m in enumerate(metabolite_ids)}
    rows, columns, values = [], [], []
    for column, rid in enumerate(ids):
        bounds.append(effective[rid])
        for metabolite, coefficient in reactions[rid]["stoichiometry"].items():
            rows.append(metabolite_index[metabolite])
            columns.append(column)
            values.append(coefficient)
    matrix = coo_matrix((values, (rows, columns)), shape=(len(metabolite_ids), len(ids))).tocsc()
    objective = np.zeros(len(ids))
    if "biomass_human" in ids:
        objective[ids.index("biomass_human")] = -1.0
    solved = linprog(
        objective,
        A_eq=matrix,
        b_eq=np.zeros(len(metabolite_ids)),
        bounds=bounds,
        method="highs",
        options={
            "time_limit": 10.0,
            "primal_feasibility_tolerance": 1e-8,
            "dual_feasibility_tolerance": 1e-8,
        },
    )
    output = {
        "success": bool(solved.success),
        "scipy_status": int(solved.status),
        "status": {0: "optimal", 1: "limit", 2: "infeasible", 3: "unbounded"}.get(
            solved.status, "error"
        ),
        "message": str(solved.message),
        "selected_reaction_count": len(ids),
    }
    if solved.success:
        output.update(
            maximum_biomass=float(-solved.fun),
            max_mass_balance_residual=float(np.max(np.abs(matrix @ solved.x), initial=0)),
        )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"preserve existing audit: {args.output}")
    reactions = read_sbml(args.human_gem)
    model_sha = sha256(args.human_gem)
    rows = []
    for study in STUDIES:
        directory = args.root / "data/processed/corneto" / study / "checkpoint_b25"
        context_path = directory / "context.json"
        context = read_json(context_path)
        context_sha = sha256(context_path)
        if context["model"]["sha256"] != model_sha:
            raise ValueError(f"model/context mismatch: {study}")
        for receipt_path in sorted((directory / "instrumented_attempts").glob("*.json")):
            receipt = read_json(receipt_path)
            condition = receipt.get("condition")
            if condition not in context["conditions"] or not receipt.get("solution"):
                continue
            if receipt.get("context_sha256") != context_sha:
                raise ValueError(f"receipt/context mismatch: {receipt_path}")
            summary = receipt["solution"]
            selected = set(summary["active_by_indicator"])
            # Positive control: retain also the reactions whose nonzero flux
            # was omitted by the indicator selection. This distinguishes a
            # bad restricted subnetwork from an SBML/LP audit implementation error.
            reported_flux_support = selected | set(summary["active_by_flux"])
            bounds = context["reaction_bounds"][condition]
            saved_growth = dict(summary["nonzero_fluxes"]).get("biomass_human", 0.0)
            rows.append(
                {
                    "study": study,
                    "condition": condition,
                    "receipt": str(receipt_path.relative_to(args.root)),
                    "receipt_sha256": sha256(receipt_path),
                    "context_sha256": context_sha,
                    "saved_biomass": saved_growth,
                    "at_context_growth_floor": restricted_lp(reactions, bounds, selected),
                    "at_saved_biomass": restricted_lp(reactions, bounds, selected, saved_growth),
                    "with_reported_flux_support_at_context_floor": restricted_lp(
                        reactions, bounds, reported_flux_support
                    ),
                }
            )
    output = {
        "schema_version": "fixed_indicator_lp_audit.v1",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "solver": "scipy.optimize.linprog/HiGHS",
        "scipy_version": scipy.__version__,
        "gurobi_sessions_requested": 0,
        "model_sha256": model_sha,
        "script_sha256": sha256(Path(__file__)),
        "results": rows,
        "claim_limit": "Restricted LP feasibility only; OCM biology and medium remain unvalidated.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp")
    temporary.write_text(json.dumps(output, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(args.output)
    print(
        json.dumps(
            {"status": "audit_completed", "attempts_audited": len(rows), "output": str(args.output)}
        )
    )


if __name__ == "__main__":
    main()
