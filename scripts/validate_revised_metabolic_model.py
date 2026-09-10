#!/usr/bin/env python3
"""Bounded GPR/pFBA revision: no b25 promotion, no calibrated-medium claim."""

import argparse
import json
from pathlib import Path

import numpy as np
from research_validation_common import load_data, sha, stamp, write_json
from run_metabolic_information_pilot import (
    ContinuousModel,
    cap_policy,
    normalize_unique_gene_ids,
    read_gprs,
    read_sbml,
)
from scipy.optimize import linprog
from scipy.sparse import bmat, eye


def validate_support(reactions, caps, fraction=0.9, seconds=30):
    model = ContinuousModel(reactions)
    maximum = model.solve(caps, set(), seconds)
    out = {"growth_maximum": maximum, "fraction": fraction}
    if maximum["status"] != "optimal" or not maximum["numerical_pass"]:
        return dict(out, status="incomplete")
    growth = maximum["maximum_biomass"]
    if growth <= 1e-8:
        return dict(out, status="no_positive_growth")
    bounds = [
        (
            max(reactions[r]["lower"], caps.get(r, (-np.inf, np.inf))[0]),
            min(reactions[r]["upper"], caps.get(r, (-np.inf, np.inf))[1]),
        )
        for r in model.ids
    ]
    n = len(bounds)
    biomass = model.ids.index("biomass_human")
    bounds[biomass] = (max(bounds[biomass][0], growth * fraction), bounds[biomass][1])
    ident = eye(n, format="csc")
    sol = linprog(
        np.r_[np.zeros(n), np.ones(n)],
        A_ub=bmat([[ident, -ident], [-ident, -ident]], format="csc"),
        b_ub=np.zeros(2 * n),
        A_eq=bmat([[model.s, model.s * 0]], format="csc"),
        b_eq=np.zeros(model.s.shape[0]),
        bounds=bounds + [(0, None)] * n,
        method="highs",
        options={
            "time_limit": seconds,
            "primal_feasibility_tolerance": 1e-8,
            "dual_feasibility_tolerance": 1e-8,
        },
    )
    out["pfba_solver_status"] = int(sol.status)
    if not sol.success:
        return dict(out, status="incomplete", message=sol.message)
    flux = sol.x[:n]
    mass = float(np.max(np.abs(model.s @ flux), initial=0))
    violation = float(
        max([0.0] + [max(lo - v, v - hi) for (lo, hi), v in zip(bounds, flux, strict=True)])
    )
    support = [r for r, v in zip(model.ids, flux, strict=True) if abs(v) > 1e-7]
    fixed = dict(caps)
    fixed.update({r: (0.0, 0.0) for r in model.ids if r not in support})
    check = model.solve(fixed, set(), seconds)
    ok = (
        mass <= 1e-6
        and violation <= 1e-6
        and check["status"] == "optimal"
        and check["numerical_pass"]
        and check["maximum_biomass"] >= growth * fraction - 1e-6
    )
    return dict(
        out,
        status="validated_support" if ok else "incomplete",
        mass_residual=mass,
        bound_violation=violation,
        pfba_biomass=float(flux[biomass]),
        total_absolute_flux=float(np.abs(flux).sum()),
        flux={r: float(v) for r, v in zip(model.ids, flux, strict=True) if abs(v) > 1e-7},
        selected_reactions=support,
        fixed_support_check=check,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--rna-root", type=Path, required=True)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    receipt = args.output / "receipt.json"
    report = {
        "status": "running",
        "started_at": stamp(),
        "cases": [],
        "script_sha256": sha(__file__),
        "model_sha256": sha(args.model),
        "claim_limit": (
            "GPR-constrained pFBA diagnostic, not minimum-cardinality CORNETO, "
            "calibrated OCMI or canonical biology."
        ),
    }
    write_json(receipt, report)
    try:
        reactions = read_sbml(args.model)
        gprs, _ = read_gprs(args.model)
        genes, _, x, metadata, sources = load_data(args.root, args.rna_root)
        genes = normalize_unique_gene_ids(genes)
        report["sources"] = sources
        source = (
            args.root
            / "data/processed/post_audit_validation_20260908/metabolic_smoke_r3/receipt.json"
        )
        prior = json.loads(source.read_text())
        assert prior["status"] == "completed" and len(prior["sample_information"]) == 4
        assert prior["model_sha256"] == report["model_sha256"]
        report["panel_source"] = {"path": str(source), "sha256": sha(source)}
        report["helper_sha256"] = {
            f: sha(Path(__file__).with_name(f))
            for f in ["run_metabolic_information_pilot.py", "research_validation_common.py"]
        }
        uptake = {r for r, v in reactions.items() if len(v["stoichiometry"]) == 1}
        report["boundary_assumption"] = (
            "Legacy model bounds; all one-metabolite boundary uptake closed "
            "only as negative control; not OCMI."
        )
        for sample in prior["sample_information"]:
            index = next(
                i for i, m in enumerate(metadata) if m["run_accession"] == sample["condition"]
            )
            expression = dict(zip(genes, x[:, index].tolist(), strict=True))
            caps, audit = cap_policy(reactions, gprs, expression, 1.0)
            negative = ContinuousModel(reactions).solve(caps, uptake, 30)
            report.setdefault("negative_controls", []).append({**sample, "result": negative})
            for fraction in [0.5, 0.9]:
                result = validate_support(reactions, caps, fraction)
                report["cases"].append({**sample, "cap_audit": audit, **result})
                write_json(receipt, report)
        passed = all(c["status"] == "validated_support" for c in report["cases"])
        negative_ok = all(
            c["result"]["status"] == "infeasible"
            or (
                c["result"]["status"] == "optimal"
                and c["result"]["numerical_pass"]
                and c["result"]["maximum_biomass"] <= 1e-8
            )
            for c in report["negative_controls"]
        )
        report["status"] = "completed" if passed and negative_ok else "incomplete"
    except Exception as exc:
        report["status"], report["error"] = "incomplete", repr(exc)
        raise
    finally:
        report["finished_at"] = stamp()
        write_json(receipt, report)
    if report["status"] != "completed":
        raise RuntimeError("revised validation gate failed; inspect receipt")


if __name__ == "__main__":
    main()
