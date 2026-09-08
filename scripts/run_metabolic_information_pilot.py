#!/usr/bin/env python3
"""Bounded, license-free LP tests of expression information and energy inputs.

The energy-uptake exclusion is a counterfactual diagnostic, NOT an OCMI
reconstruction. GPR-derived caps are uncalibrated sensitivity assumptions.
This script never changes frozen b25 contexts or promotes a canonical receipt.
"""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import numpy as np
from audit_metabolic_scientific_validity import CORE, FBC, read_sbml
from research_validation_common import STUDIES, load_data, sha, stamp, write_json
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def decode(identifier, prefix):
    return re.sub(r"__(\d+)__", lambda m: chr(int(m[1])), identifier.removeprefix(prefix))


def normalize_gene_id(identifier):
    """Remove an Ensembl gene version without collapsing a PAR_Y identity."""
    return re.sub(r"^(ENSG\d+)\.\d+(_PAR_Y)?$", r"\1\2", identifier)


def normalize_unique_gene_ids(identifiers):
    normalized = [normalize_gene_id(identifier) for identifier in identifiers]
    duplicates = sorted(k for k, n in Counter(normalized).items() if n > 1)
    if duplicates:
        raise ValueError("gene ID normalization creates duplicates: " + ", ".join(duplicates[:10]))
    return normalized


def read_gprs(path):
    model = ET.parse(path).getroot().find(CORE + "model")
    products = model.find(FBC + "listOfGeneProducts")
    gene_labels = {
        p.attrib[FBC + "id"]: p.attrib.get(FBC + "label", p.attrib[FBC + "id"]) for p in products
    }
    gene_labels = {k: normalize_gene_id(decode(v, "G_")) for k, v in gene_labels.items()}

    def tree(element):
        tag = element.tag.split("}")[-1]
        if tag == "geneProductRef":
            return ("gene", gene_labels[element.attrib[FBC + "geneProduct"]])
        if tag not in ("and", "or", "geneProductAssociation"):
            raise ValueError("unknown GPR operation: " + tag)
        children = [tree(c) for c in element]
        if tag == "geneProductAssociation":
            if len(children) != 1:
                raise ValueError("invalid GPR root")
            return children[0]
        return (tag, children)

    result = {}
    for reaction in model.find(CORE + "listOfReactions"):
        gpr = reaction.find(FBC + "geneProductAssociation")
        if gpr is not None:
            result[decode(reaction.attrib["id"], "R_")] = tree(gpr)
    return result, sorted(set(gene_labels.values()))


def evaluate(rule, expression):
    if rule[0] == "gene":
        return expression.get(rule[1])
    values = [evaluate(child, expression) for child in rule[1]]
    # Any missing subunit/isoenzyme makes capacity unknown, never zero.
    if not values or any(v is None for v in values):
        return None
    return min(values) if rule[0] == "and" else sum(values)


def cap_policy(reactions, gprs, expression, scale):
    caps, audit = {}, Counter()
    for rid, rule in gprs.items():
        value = evaluate(rule, expression)
        if value is None:
            audit["missing_gene_capacity_unknown"] += 1
            continue
        if value == 0:
            audit["zero_rna_not_forced_to_knockout"] += 1
            continue
        if not np.isfinite(value) or value < 0:
            raise ValueError("invalid GPR expression")
        capacity = float(value * scale)
        reaction = reactions[rid]
        # Bounds preserve reverse-only and irreversible reaction directions.
        caps[rid] = (max(reaction["lower"], -capacity), min(reaction["upper"], capacity))
    audit["positive_caps"] = len(caps)
    return caps, dict(audit)


def block_uptake(reaction):
    stoich = reaction["stoichiometry"]
    if len(stoich) != 1:
        raise ValueError("uptake direction requires a one-metabolite exchange")
    coefficient = next(iter(stoich.values()))
    if coefficient == 0:
        raise ValueError("zero exchange coefficient")
    lower, upper = reaction["lower"], reaction["upper"]
    return (max(lower, 0.0), upper) if coefficient < 0 else (lower, min(upper, 0.0))


class ContinuousModel:
    def __init__(self, reactions):
        self.reactions = reactions
        self.ids = sorted(reactions)
        mids = sorted({m for r in reactions.values() for m in r["stoichiometry"]})
        lookup = {m: i for i, m in enumerate(mids)}
        rows, columns, values = [], [], []
        for j, rid in enumerate(self.ids):
            for mid, coefficient in reactions[rid]["stoichiometry"].items():
                rows.append(lookup[mid])
                columns.append(j)
                values.append(coefficient)
        self.s = coo_matrix((values, (rows, columns)), shape=(len(mids), len(self.ids))).tocsc()
        self.c = np.zeros(len(self.ids))
        self.c[self.ids.index("biomass_human")] = -1.0

    def solve(self, caps, blocked, seconds):
        bounds = []
        for rid in self.ids:
            r = self.reactions[rid]
            lower, upper = r["lower"], r["upper"]
            if rid in caps:
                lower, upper = max(lower, caps[rid][0]), min(upper, caps[rid][1])
            if rid in blocked:
                lo, hi = block_uptake(r)
                lower, upper = max(lower, lo), min(upper, hi)
            if lower > upper:
                return {"status": "infeasible_bounds", "reaction": rid}
            bounds.append((lower, upper))
        solution = linprog(
            self.c,
            A_eq=self.s,
            b_eq=np.zeros(self.s.shape[0]),
            bounds=bounds,
            method="highs",
            options={
                "time_limit": seconds,
                "primal_feasibility_tolerance": 1e-8,
                "dual_feasibility_tolerance": 1e-8,
            },
        )
        result = {
            "status": {0: "optimal", 1: "limit", 2: "infeasible", 3: "unbounded"}.get(
                solution.status, "error"
            ),
            "solver_message": solution.message,
        }
        if solution.success:
            flux = dict(zip(self.ids, solution.x.tolist(), strict=True))
            mass_error = float(np.max(np.abs(self.s @ solution.x), initial=0))
            bound_error = max(
                [0.0]
                + [max(lo - v, v - hi) for (lo, hi), v in zip(bounds, solution.x, strict=True)]
            )
            result.update(
                maximum_biomass=float(-solution.fun),
                max_mass_balance_residual=mass_error,
                max_bound_violation=float(bound_error),
                numerical_pass=mass_error <= 1e-6 and bound_error <= 1e-6,
                flux_support_count=int(np.sum(np.abs(solution.x) > 1e-7)),
                expression_caps_with_nonzero_flux=sum(abs(flux[r]) > 1e-7 for r in caps),
                diagnostic_exchanges={
                    r: flux.get(r)
                    for r in ["EX_atp[e]", "EX_pep[e]", "EX_pcreat[e]", "EX_glc_D[e]"]
                },
            )
        return result


def choose_conditions(context, count):
    # Response-blind coverage of lowest and median cap counts, not outcomes.
    candidates = sorted(
        context["conditions"], key=lambda c: (len(context["reaction_bounds"][c]), c)
    )
    if count == 1:
        return [candidates[len(candidates) // 2]]
    if count == 2:
        return [candidates[0], candidates[len(candidates) // 2]]
    raise ValueError("bounded pilot supports only 1 or 2 OCM per study")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--rna-root", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--per-study", type=int, default=2)
    parser.add_argument("--shuffles", type=int, default=10)
    parser.add_argument("--lp-seconds", type=int, default=10)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "metabolic_information_pilot.v1",
        "status": "running",
        "started_at_utc": stamp(),
        "script_sha256": sha(__file__),
        "common_code_sha256": sha(Path(__file__).with_name("research_validation_common.py")),
        "model_sha256": sha(args.human_gem),
        "gurobi_sessions_requested": 0,
        "results": [],
        "settings": {
            "per_study": args.per_study,
            "shuffles": args.shuffles,
            "lp_seconds": args.lp_seconds,
            "real_gpr_scales": [0.1, 1.0, 10.0],
            "shuffle_scale": 1.0,
            "seed": 20260908,
        },
        "claim_limit": (
            "LP input-sensitivity diagnostic only. No calibrated OCMI medium, measured "
            "flux/capacity, sparse selection, canonical receipt, or TPI1 interpretation. "
            "Old growth floor is not imposed on the alternative scenarios."
        ),
    }
    output = args.output_dir / "receipt.json"
    write_json(output, report)
    try:
        reactions = read_sbml(args.human_gem)
        gprs, model_genes = read_gprs(args.human_gem)
        genes, names, x, metadata, sources = load_data(args.root, args.rna_root)
        normalized_ids = normalize_unique_gene_ids(genes)
        report["sources"] = sources
        run_index = {m["run_accession"]: i for i, m in enumerate(metadata)}
        model = ContinuousModel(reactions)
        energy = {"EX_atp[e]", "EX_pep[e]", "EX_pcreat[e]"}
        if not energy.issubset(reactions):
            raise ValueError("diagnostic exchange identifiers missing")
        for study in STUDIES:
            path = args.root / "data/processed/corneto" / study / "checkpoint_b25/context.json"
            context = json.loads(path.read_text())
            if context["model"]["sha256"] != report["model_sha256"]:
                raise ValueError("frozen model hash mismatch")
            report.setdefault("contexts", []).append({"path": str(path), "sha256": sha(path)})
            for condition in choose_conditions(context, args.per_study):
                expression = dict(
                    zip(normalized_ids, x[:, run_index[condition]].tolist(), strict=True)
                )
                covered = sorted(set(model_genes) & set(expression))
                if len(covered) < 0.8 * len(model_genes):
                    raise ValueError("less than 80% model-gene identity coverage; inspect mapping")
                report.setdefault("sample_information", []).append(
                    {
                        "study": study,
                        "condition": condition,
                        "model_genes_covered": len(covered),
                        "model_gene_count": len(model_genes),
                    }
                )
                old_caps = {
                    r: b
                    for r, b in context["reaction_bounds"][condition].items()
                    if r != "biomass_human"
                }
                policies = [("expression_null", {}, {}), ("frozen_b25", old_caps, {})]
                for scale in (0.1, 1.0, 10.0):
                    caps, audit = cap_policy(reactions, gprs, expression, scale)
                    policies.append((f"gpr_real_scale_{scale}", caps, audit))
                rng = np.random.default_rng(20260908)
                for iteration in range(args.shuffles):
                    # Permute values among measured MODEL genes; preserve zero counts
                    # and the expression marginal, but destroy gene/value identity.
                    shuffled = dict(expression)
                    values = rng.permutation([expression[g] for g in covered])
                    shuffled.update(zip(covered, values.tolist(), strict=True))
                    caps, audit = cap_policy(reactions, gprs, shuffled, 1.0)
                    policies.append((f"gpr_shuffled_{iteration}", caps, audit))
                for label, caps, audit in policies:
                    for medium, blocked in (
                        ("legacy_default", set()),
                        ("exclude_three_energy_uptakes", energy),
                    ):
                        result = model.solve(caps, blocked, args.lp_seconds)
                        report["results"].append(
                            {
                                "study": study,
                                "condition": condition,
                                "policy": label,
                                "medium_diagnostic": medium,
                                "expression_cap_count": len(caps),
                                "gpr_audit": audit,
                                **result,
                            }
                        )
                        write_json(output, report)
                print(study, condition, "LP diagnostic recorded", flush=True)
        errors = [
            r
            for r in report["results"]
            if r["status"] not in ("optimal", "infeasible", "infeasible_bounds")
            or (r["status"] == "optimal" and not r["numerical_pass"])
        ]
        report["status"] = "completed" if not errors else "incomplete"
        report["unresolved_solver_case_count"] = len(errors)
    except Exception as exc:
        report["status"], report["error"] = "incomplete", repr(exc)
        raise
    finally:
        report["finished_at_utc"] = stamp()
        write_json(output, report)
    if report["status"] != "completed":
        raise RuntimeError("pilot has unresolved numerical/solver cases")
    print("completed", output, flush=True)


if __name__ == "__main__":
    main()
