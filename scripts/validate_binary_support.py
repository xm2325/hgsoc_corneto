"""Small topology-restricted binary linking test; not full CORNETO optimization."""

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
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import bmat, diags, lil_matrix


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--rna", type=Path, required=True)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--disable-presolve", action="store_true")
    p.add_argument("--native-indicators", action="store_true")
    p.add_argument("--case-index", type=int, help="Explicit index in the audited pFBA panel")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    report = {
        "status": "running",
        "started_at": stamp(),
        "script_sha256": sha(__file__),
        "claim_limit": (
            "One pFBA-support-restricted binary linking test, "
            "not global CORNETO optimum or biology."
        ),
    }
    target = a.output / "receipt.json"
    write_json(target, report)
    try:
        source = a.root / "data/processed/revised_gpr_pfba_20260910/receipt.json"
        prior = json.loads(source.read_text())
        assert prior["status"] == "completed" and prior["model_sha256"] == sha(a.model)
        case_index = a.case_index
        if case_index is None:
            case_index = next(i for i, c in enumerate(prior["cases"]) if c["fraction"] == 0.9)
        if not 0 <= case_index < len(prior["cases"]):
            raise ValueError("case index outside the audited pFBA panel")
        case = prior["cases"][case_index]
        report.update(
            case_index=case_index,
            growth_fraction=case["fraction"],
            source_sha256=sha(source),
            model_sha256=sha(a.model),
            condition=case["condition"],
            maximum_seconds=120,
            requested_gap=1e-4,
        )
        reactions = read_sbml(a.model)
        gprs, _ = read_gprs(a.model)
        genes, _, x, meta, sources = load_data(a.root, a.rna)
        report["sources"] = sources
        index = next(i for i, m in enumerate(meta) if m["run_accession"] == case["condition"])
        expr = dict(zip(normalize_unique_gene_ids(genes), x[:, index].tolist(), strict=True))
        caps, _ = cap_policy(reactions, gprs, expr, 1.0)
        subset = {r: reactions[r] for r in case["selected_reactions"]}
        model = ContinuousModel(subset)
        n = len(model.ids)
        lo = np.array(
            [max(subset[r]["lower"], caps.get(r, (-np.inf, np.inf))[0]) for r in model.ids]
        )
        hi = np.array(
            [min(subset[r]["upper"], caps.get(r, (-np.inf, np.inf))[1]) for r in model.ids]
        )
        growth = case["growth_maximum"]["maximum_biomass"] * case["fraction"]
        row = lil_matrix((1, 2 * n))
        row[0, model.ids.index("biomass_human")] = 1
        ident = diags(np.ones(n))
        constraints = [
            LinearConstraint(bmat([[model.s, model.s * 0]]), 0, 0),
            LinearConstraint(bmat([[ident, -diags(hi)]]), -np.inf, 0),
            LinearConstraint(bmat([[ident, -diags(lo)]]), 0, np.inf),
            LinearConstraint(row.tocsc(), growth, np.inf),
        ]
        witness = np.r_[[case["flux"][r] for r in model.ids], np.ones(n)]
        violations = []
        for constraint in constraints:
            activity = constraint.A @ witness
            violations.append(
                float(max(0, np.max(constraint.lb - activity), np.max(activity - constraint.ub)))
            )
        lower = np.r_[np.minimum(lo, 0), np.zeros(n)]
        upper = np.r_[np.maximum(hi, 0), np.ones(n)]
        witness_bounds = float(max(0, np.max(lower - witness), np.max(witness - upper)))
        report["all_on_witness"] = {
            "constraint_violations": violations,
            "bound_violation": witness_bounds,
        }
        report["presolve"] = not a.disable_presolve
        write_json(target, report)
        if max(violations + [witness_bounds]) > 1e-6:
            raise ValueError("known flux is not feasible in the assembled MILP")
        if a.native_indicators:
            from native_indicator_solver import solve_native

            report["native_helper_sha256"] = sha(
                Path(__file__).with_name("native_indicator_solver.py")
            )
            report["solver_policy"] = {
                "engine": "gurobi_native_indicators",
                "FeasibilityTol": 1e-9,
                "IntFeasTol": 1e-9,
                "IntegralityFocus": 1,
                "warm_start": "saved_pFBA_flux_and_all_indicators_one",
            }
            write_json(target, report)
            sol = solve_native(
                model.s, lo, hi, model.ids.index("biomass_human"), growth, witness, a.output
            )
        else:
            sol = milp(
                np.r_[np.zeros(n), np.ones(n)],
                integrality=np.r_[np.zeros(n), np.ones(n)],
                bounds=Bounds(
                    np.r_[np.minimum(lo, 0), np.zeros(n)], np.r_[np.maximum(hi, 0), np.ones(n)]
                ),
                constraints=constraints,
                options={
                    "time_limit": 120,
                    "mip_rel_gap": 1e-4,
                    "disp": True,
                    "presolve": not a.disable_presolve,
                },
            )
        report.update(solver_status=int(sol.status), message=sol.message, restricted_reactions=n)
        for key in ["fun", "mip_gap", "mip_dual_bound", "mip_node_count"]:
            v = getattr(sol, key, None)
            report[key] = float(v) if v is not None and np.isfinite(v) else None
        report["status"] = "no_incumbent"
        if sol.x is not None:
            v, y = sol.x[:n], sol.x[n:]
            selected = {r for r, b in zip(model.ids, y, strict=True) if b > 0.5}
            fixed = dict(caps)
            fixed.update({r: (0, 0) for r in reactions if r not in selected})
            check = ContinuousModel(reactions).solve(fixed, set(), 30)
            residual = float(np.max(np.abs(model.s @ v), initial=0))
            integer = float(np.max(np.abs(y - np.round(y)), initial=0))
            leakage = max(
                [0.0]
                + [abs(float(f)) for r, f in zip(model.ids, v, strict=True) if r not in selected]
            )
            link_violation = float(max(0, np.max(v - hi * y), np.max(lo * y - v)))
            growth_violation = max(0.0, growth - float(v[model.ids.index("biomass_human")]))
            ok = (
                check["status"] == "optimal"
                and check["numerical_pass"]
                and check["maximum_biomass"] >= growth - 1e-6
                and residual <= 1e-6
                and integer <= 1e-6
                and leakage <= 1e-7
                and link_violation <= 1e-6
                and growth_violation <= 1e-6
            )
            report.update(
                selected_reactions=sorted(selected),
                fixed_support_check=check,
                mass_residual=residual,
                integrality_error=integer,
                unselected_flux=leakage,
                link_violation=link_violation,
                growth_violation=growth_violation,
                flux={r: float(f) for r, f in zip(model.ids, v, strict=True)},
                status="validated_restricted_optimum"
                if ok and sol.status == 0
                else "validated_partial_incumbent"
                if ok
                else "invalid_incumbent",
            )
    except Exception as exc:
        report["status"], report["error"] = "incomplete", repr(exc)
        raise
    finally:
        report["finished_at"] = stamp()
        write_json(target, report)
    if report["status"] != "validated_restricted_optimum":
        raise RuntimeError("binary scientific gate not passed; receipt preserved")


if __name__ == "__main__":
    main()
