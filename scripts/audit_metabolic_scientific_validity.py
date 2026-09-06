#!/usr/bin/env python3
"""Read-only model/input/receipt audit; no optimization or licence is required.

Checks sparse saved summaries, not an entire MIP search tree. Truncated flux
summaries can omit values below their reporting tolerance; residuals are
therefore diagnostics, not a full-precision feasibility certificate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import re
import statistics
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

STUDIES = ("E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568")
CORE = "{http://www.sbml.org/sbml/level3/version1/core}"
FBC = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def read_sbml(path):
    model = ET.parse(path).getroot().find(f"{CORE}model")
    params = {
        p.attrib["id"]: float(p.attrib["value"]) for p in model.find(f"{CORE}listOfParameters")
    }
    reactions = {}
    for reaction in model.find(f"{CORE}listOfReactions"):
        rid = reaction.attrib["id"]
        rid = re.sub(r"__(\d+)__", lambda m: chr(int(m[1])), rid.removeprefix("R_"))
        stoich = {}
        for group, sign in (("listOfReactants", -1), ("listOfProducts", 1)):
            refs = reaction.find(f"{CORE}{group}")
            for ref in refs if refs is not None else []:
                sid = ref.attrib["species"]
                stoich[sid] = stoich.get(sid, 0.0) + sign * float(ref.attrib["stoichiometry"])
        reactions[rid] = {
            "lower": params[reaction.attrib[f"{FBC}lowerFluxBound"]],
            "upper": params[reaction.attrib[f"{FBC}upperFluxBound"]],
            "stoichiometry": stoich,
        }
    return reactions


def range_summary(values):
    values = list(values)
    return (
        {"min": min(values), "median": statistics.median(values), "max": max(values)}
        if values
        else None
    )


def bound_violations(flux, bounds, tolerance=1e-6):
    return [
        rid
        for rid, (lower, upper) in bounds.items()
        if flux.get(rid, 0.0) < lower - tolerance or flux.get(rid, 0.0) > upper + tolerance
    ]


def inspect_solution(solution, reactions, sample_bounds, objective):
    flux = dict(solution["nonzero_fluxes"])
    if any(not math.isfinite(float(x)) for x in flux.values()):
        raise ValueError("nonfinite saved flux")
    unknown = sorted(set(flux) - set(reactions))
    residuals = {}
    for rid, value in flux.items():
        for sid, coefficient in reactions.get(rid, {}).get("stoichiometry", {}).items():
            residuals[sid] = residuals.get(sid, 0.0) + value * coefficient
    indicators = set(solution["active_by_indicator"])
    mismatch = sorted(set(flux) - indicators)
    caps = set(sample_bounds) - set(objective)
    used_caps = sorted(r for r in caps if abs(flux.get(r, 0.0)) > 1e-7)
    model_bounds = {r: (v["lower"], v["upper"]) for r, v in reactions.items()}
    return {
        "unknown_reactions": unknown,
        "max_sparse_summary_mass_balance_residual": max(map(abs, residuals.values()), default=0.0),
        "model_bound_violation_count": len(bound_violations(flux, model_bounds)),
        "context_bound_violation_count": len(bound_violations(flux, sample_bounds)),
        "flux_without_selected_indicator_count": len(mismatch),
        "flux_without_selected_indicator_magnitude": range_summary(abs(flux[r]) for r in mismatch),
        "selected_indicator_count": len(indicators),
        "nonzero_reported_flux_count": len(flux),
        "expression_cap_count": len(caps),
        "expression_caps_with_reported_flux": used_caps,
        "biomass_flux": flux.get("biomass_human", 0.0),
        "diagnostic_exchange_fluxes": {
            r: flux.get(r, 0.0) for r in ("EX_atp[e]", "EX_pep[e]", "EX_pcreat[e]", "EX_glc_D[e]")
        },
    }


def snapshot_sources(root, snapshot=None):
    snapshot = snapshot or root / "evidence/roihu_result_snapshot.json"
    results = []

    def visit(value):
        if isinstance(value, dict):
            if isinstance(value.get("path"), str) and isinstance(value.get("sha256"), str):
                source = root / value["path"]
                actual = sha256(source) if source.is_file() else None
                results.append(
                    {
                        "path": value["path"],
                        "expected_sha256": value["sha256"],
                        "actual_sha256": actual,
                        "matches": actual == value["sha256"],
                    }
                )
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    if snapshot.is_file():
        visit(read_json(snapshot))
    return results


def audit(root, model_path, registry=None, snapshot=None):
    reactions = read_sbml(model_path)
    model_sha = sha256(model_path)
    all_bounds, contexts, cohort_rows, attempts = {}, {}, [], []
    for study in STUDIES:
        checkpoint = root / "data/processed/corneto" / study / "checkpoint_b25"
        context_path = checkpoint / "context.json"
        context = read_json(context_path)
        context_sha = sha256(context_path)
        if context["model"]["sha256"] != model_sha:
            raise ValueError(f"model hash mismatch: {study}")
        if len(context["conditions"]) != len(set(context["conditions"])):
            raise ValueError(f"duplicate conditions: {study}")
        contexts[study] = context
        caps = []
        for condition in context["conditions"]:
            if condition in all_bounds:
                raise ValueError(f"duplicate condition across contexts: {condition}")
            all_bounds[condition] = context["reaction_bounds"][condition]
            caps.append(len(set(all_bounds[condition]) - set(context["objectives"][condition])))
        canonical_files = sorted((checkpoint / "independent").glob("*.json"))
        cohort_rows.append(
            {
                "study": study,
                "context_path": str(context_path.relative_to(root)),
                "context_sha256": context_sha,
                "condition_count": len(context["conditions"]),
                "canonical_file_count_unvalidated": len(canonical_files),
                "expression_caps_per_condition": range_summary(caps),
                "conditions_with_zero_expression_caps": sum(v == 0 for v in caps),
                "distinct_bound_contexts": len(
                    {json.dumps(all_bounds[c], sort_keys=True) for c in context["conditions"]}
                ),
                "growth_optimum": range_summary(context["objective"]["growth_optima"].values()),
                "objective": {
                    k: context["objective"][k]
                    for k in ("independent_lambda", "joint_lambda", "growth_fraction")
                },
                "candidate_reaction_ids": context["candidate_selection"]["selected_reaction_ids"],
            }
        )
        for path in sorted((checkpoint / "instrumented_attempts").glob("*.json")):
            receipt = read_json(path)
            condition = receipt.get("condition")
            if condition not in context["conditions"]:
                continue  # joint attempts are a different evidence object
            instrumentation = receipt.get("instrumentation", {})
            artifacts = instrumentation.get("artifacts", {})
            artifact_checks = {}
            for label in ("solution", "mip_start", "gurobi_log"):
                p = Path(artifacts[label]) if artifacts.get(label) else None
                artifact_checks[label] = {
                    "exists_nonempty": bool(p and p.is_file() and p.stat().st_size),
                    "sha256": sha256(p) if p and p.is_file() else None,
                }
            row = {
                "study": study,
                "condition": condition,
                "array_index": receipt.get("array_index"),
                "receipt_path": str(path.relative_to(root)),
                "receipt_sha256": sha256(path),
                "status": receipt.get("status"),
                "scientific_success": receipt.get("scientific_success"),
                "context_matches": receipt.get("context_sha256") == context_sha,
                "telemetry": instrumentation.get("telemetry"),
                "artifact_checks": artifact_checks,
            }
            if receipt.get("solution") and row["context_matches"]:
                row["diagnostics"] = inspect_solution(
                    receipt["solution"],
                    reactions,
                    all_bounds[condition],
                    context["objectives"][condition],
                )
                row["_flux"] = dict(receipt["solution"]["nonzero_fluxes"])
            attempts.append(row)
    for row in attempts:
        if "_flux" in row:
            flux = row.pop("_flux")
            feasible = [c for c, bounds in all_bounds.items() if not bound_violations(flux, bounds)]
            row["diagnostics"]["contexts_accepting_saved_flux_bounds"] = len(feasible)
            row["diagnostics"]["contexts_tested"] = len(all_bounds)
    overlaps = []
    for a, b in itertools.combinations(cohort_rows, 2):
        x, y = set(a["candidate_reaction_ids"]), set(b["candidate_reaction_ids"])
        overlaps.append(
            {
                "studies": [a["study"], b["study"]],
                "intersection": len(x & y),
                "union": len(x | y),
                "jaccard": len(x & y) / len(x | y),
            }
        )
    registry = registry or root / "evidence/study_ocm_registry.tsv"
    registry_result = None
    if registry.is_file():
        with registry.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        primary = [r for r in rows if r["primary_cohort_eligible"].lower() == "true"]
        families = Counter(r["patient_id"] for r in primary)
        registry_result = {
            "sha256": sha256(registry),
            "all_rows": len(rows),
            "primary_rows": len(primary),
            "patients": len(families),
            "primary_by_study": dict(Counter(r["study_accession"] for r in primary)),
            "repeated_patients": {k: v for k, v in families.items() if v > 1},
        }
    regulatory = []
    grid = root / "data/processed/corneto/regulatory_multisample/v1/grid"
    for path in sorted(grid.glob("*/*.json")):
        r = read_json(path)
        regulatory.append(
            {
                "path": str(path.relative_to(root)),
                "sha256": sha256(path),
                "status": r.get("status"),
                "solver": r.get("solver"),
                "processed_graph": r.get("processed_graph"),
                "scope_counts": r.get("scope_counts"),
                "union_count": r.get("selected_edge_union_count"),
                "lambda": r.get("method", {}).get("lambda_nominal"),
            }
        )
    return {
        "schema_version": "hgsoc_scientific_validity_audit.v1",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "model_sha256": model_sha,
        "model_reaction_count": len(reactions),
        "sample_registry": registry_result,
        "cohorts": cohort_rows,
        "candidate_overlaps": overlaps,
        "attempts": attempts,
        "regulatory_grid": regulatory,
        "historical_snapshot_source_verification": snapshot_sources(root, snapshot),
        "interpretation": "Diagnostics only; no new solver result or biological validation.",
        "limitations": [
            "Saved flux summaries omit values below reporting tolerance.",
            "Cross-context bound compatibility does not prove identical feasible sets.",
            "No fixed-indicator LP, medium calibration or null-model test was solved.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--snapshot", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            audit(args.root, args.human_gem, args.registry, args.snapshot),
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
