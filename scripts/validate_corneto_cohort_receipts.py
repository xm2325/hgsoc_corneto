#!/usr/bin/env python3
"""Validate and compare response-blind CORNETO cohort receipts.

This is a read-only post-processing gate.  It never loads SBML, CORNETO,
CVXPY, or a solver license.  Each ``--receipt`` argument has the form
``STUDY=PATH``.  A malformed or incomplete receipt makes the process fail
closed and no output file is written.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any


class ReceiptError(ValueError):
    """Raised when a receipt is not structurally usable."""


EXPECTED_SAMPLE_COUNTS = {
    "E-MTAB-7223": 9,
    "E-MTAB-10801": 13,
    "E-MTAB-11000": 11,
    "E-MTAB-14568": 27,
}
EXPECTED_SOLVER = "gurobi"
EXPECTED_CANDIDATES = 25
EXPECTED_GROWTH_FRACTION = 0.9
EXPECTED_INDEPENDENT_LAMBDA = 0.1
EXPECTED_JOINT_LAMBDA = 1.0
ACCEPTED_STATUSES = {"optimal", "optimal_inaccurate"}


def _obj(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ReceiptError(f"{name} must be an object")
    return value


def _list(value: Any, name: str) -> list[Any]:
    if not isinstance(value, list):
        raise ReceiptError(f"{name} must be an array")
    return value


def _strings(value: Any, name: str) -> list[str]:
    values = _list(value, name)
    if any(not isinstance(item, str) or not item for item in values):
        raise ReceiptError(f"{name} must contain non-empty strings")
    if len(values) != len(set(values)):
        raise ReceiptError(f"{name} contains duplicates")
    return list(values)


def _number(value: Any, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool):
        raise ReceiptError(f"{name} must be numeric")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ReceiptError(f"{name} must be numeric") from error
    if not math.isfinite(parsed):
        raise ReceiptError(f"{name} must be finite")
    if positive and parsed <= 0:
        raise ReceiptError(f"{name} must be positive")
    return parsed


def _integer(value: Any, name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ReceiptError(f"{name} must be an integer >= {minimum}")
    return value


def _parse_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise ReceiptError(f"receipt must be STUDY=PATH, got {spec!r}")
    study, raw_path = spec.split("=", 1)
    if not study or not raw_path:
        raise ReceiptError(f"receipt must be STUDY=PATH, got {spec!r}")
    return study, Path(raw_path)


def _jaccard(left: set[str], right: set[str]) -> float | None:
    union = left | right
    return None if not union else len(left & right) / len(union)


def _read_receipt(study: str, path: Path) -> dict[str, Any]:
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReceiptError(f"{study}: cannot read JSON: {error}") from error
    root = _obj(receipt, f"{study} receipt")
    if root.get("status") != "completed":
        raise ReceiptError(f"{study}: status={root.get('status')!r}, expected 'completed'")
    if root.get("study_accession") != study:
        raise ReceiptError(
            f"{study}: study_accession={root.get('study_accession')!r} does not match"
        )
    if study not in EXPECTED_SAMPLE_COUNTS:
        raise ReceiptError(f"{study}: study is not in the frozen four-cohort contract")
    if root.get("primary_only") is not True:
        raise ReceiptError(f"{study}: primary_only must be true")
    if not isinstance(root.get("repo_commit"), str) or not root["repo_commit"]:
        raise ReceiptError(f"{study}: repo_commit is missing")
    samples = _list(root.get("samples"), f"{study}.samples")
    sample_ids: list[str] = []
    for index, sample in enumerate(samples):
        item = _obj(sample, f"{study}.samples[{index}]")
        run_id = item.get("run_accession")
        if not isinstance(run_id, str) or not run_id:
            raise ReceiptError(f"{study}.samples[{index}].run_accession is invalid")
        sample_ids.append(run_id)
    if len(sample_ids) != len(set(sample_ids)):
        raise ReceiptError(f"{study}: duplicate sample run_accession")
    if root.get("sample_count") != len(sample_ids):
        raise ReceiptError(f"{study}: sample_count disagrees with samples")
    if len(sample_ids) != EXPECTED_SAMPLE_COUNTS[study]:
        raise ReceiptError(
            f"{study}: sample_count={len(sample_ids)}, expected {EXPECTED_SAMPLE_COUNTS[study]}"
        )

    solver = _obj(root.get("solver"), f"{study}.solver")
    requested = solver.get("requested")
    used = solver.get("used")
    if str(requested).casefold() != EXPECTED_SOLVER or str(used).casefold() != EXPECTED_SOLVER:
        raise ReceiptError(
            f"{study}: solver requested/used must both be {EXPECTED_SOLVER!r}"
        )
    if solver.get("fallback_reason") is not None or solver.get("solve_fallback") is not None:
        raise ReceiptError(f"{study}: solver fallback is not permitted")

    candidate = _obj(root.get("candidate_selection"), f"{study}.candidate_selection")
    candidate_ids = set(_strings(candidate.get("selected_reaction_ids"), f"{study}.candidate IDs"))
    if candidate.get("selected_count") != len(candidate_ids):
        raise ReceiptError(f"{study}: selected_count disagrees with selected_reaction_ids")
    if candidate.get("max_candidates") != EXPECTED_CANDIDATES:
        raise ReceiptError(
            f"{study}: max_candidates={candidate.get('max_candidates')!r}, "
            f"expected {EXPECTED_CANDIDATES}"
        )
    if not candidate_ids or len(candidate_ids) > EXPECTED_CANDIDATES:
        raise ReceiptError(f"{study}: selected candidate count is outside 1..{EXPECTED_CANDIDATES}")

    objective = _obj(root.get("objective"), f"{study}.objective")
    expected_objective = {
        "growth_fraction": EXPECTED_GROWTH_FRACTION,
        "independent_lambda": EXPECTED_INDEPENDENT_LAMBDA,
        "joint_lambda": EXPECTED_JOINT_LAMBDA,
    }
    for field, expected in expected_objective.items():
        actual = _number(objective.get(field), f"{study}.objective.{field}")
        if not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
            raise ReceiptError(f"{study}: objective.{field}={actual}, expected {expected}")
    if _number(objective.get("biomass_coefficient"), f"{study}.objective.biomass_coefficient") != -1:
        raise ReceiptError(f"{study}: biomass_coefficient must be -1")
    optima = _obj(objective.get("growth_optima"), f"{study}.objective.growth_optima")
    targets = _obj(objective.get("growth_targets"), f"{study}.objective.growth_targets")
    if set(optima) != set(sample_ids) or set(targets) != set(sample_ids):
        raise ReceiptError(f"{study}: growth optimum/target keys disagree with samples")
    for run_id in sample_ids:
        optimum = _number(optima[run_id], f"{study}.growth_optima[{run_id}]", positive=True)
        target = _number(targets[run_id], f"{study}.growth_targets[{run_id}]", positive=True)
        if not math.isclose(
            target,
            EXPECTED_GROWTH_FRACTION * optimum,
            rel_tol=1e-9,
            abs_tol=1e-10,
        ):
            raise ReceiptError(f"{study}: growth target is inconsistent for {run_id}")

    corneto = _obj(root.get("corneto"), f"{study}.corneto")
    if str(corneto.get("solver")).casefold() != EXPECTED_SOLVER:
        raise ReceiptError(f"{study}: corneto.solver must be {EXPECTED_SOLVER!r}")
    for field, expected in (
        ("independent_lambda", EXPECTED_INDEPENDENT_LAMBDA),
        ("joint_lambda", EXPECTED_JOINT_LAMBDA),
    ):
        actual = _number(corneto.get(field), f"{study}.corneto.{field}")
        if not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
            raise ReceiptError(f"{study}: corneto.{field}={actual}, expected {expected}")
    active_tolerance = _number(
        corneto.get("active_tolerance"), f"{study}.corneto.active_tolerance", positive=True
    )
    if corneto.get("conditions") != sample_ids:
        raise ReceiptError(f"{study}: corneto.conditions disagrees with samples")
    status_counts: dict[str, dict[str, int]] = {}
    active_unions: dict[str, set[str]] = {}
    for label in ("independent", "joint"):
        summaries = _list(corneto.get(label), f"{study}.corneto.{label}")
        if len(summaries) != len(sample_ids):
            raise ReceiptError(f"{study}: {label} summary length disagrees with samples")
        union: set[str] = set()
        counts: Counter[str] = Counter()
        for index, summary in enumerate(summaries):
            item = _obj(summary, f"{study}.{label}[{index}]")
            if item.get("condition") != sample_ids[index]:
                raise ReceiptError(f"{study}.{label}[{index}].condition disagrees with samples")
            status = item.get("status")
            if not isinstance(status, str) or not status:
                raise ReceiptError(f"{study}.{label}[{index}].status is missing")
            if status.casefold() not in ACCEPTED_STATUSES:
                raise ReceiptError(f"{study}.{label}[{index}].status={status!r} is not optimal")
            _number(
                item.get("problem_objective_value"),
                f"{study}.{label}[{index}].problem_objective_value",
            )
            counts[status] += 1
            active = _strings(item.get("active_by_flux"), f"{study}.{label}[{index}].active_by_flux")
            nonzero = _list(
                item.get("nonzero_fluxes"), f"{study}.{label}[{index}].nonzero_fluxes"
            )
            nonzero_ids: set[str] = set()
            for flux_index, pair in enumerate(nonzero):
                if not isinstance(pair, list) or len(pair) != 2:
                    raise ReceiptError(
                        f"{study}.{label}[{index}].nonzero_fluxes[{flux_index}] must be [id, value]"
                    )
                reaction_id, flux = pair
                if not isinstance(reaction_id, str) or not reaction_id or reaction_id in nonzero_ids:
                    raise ReceiptError(f"{study}.{label}[{index}] has invalid/duplicate flux ID")
                nonzero_ids.add(reaction_id)
                if abs(
                    _number(flux, f"{study}.{label}[{index}].nonzero_fluxes[{flux_index}][1]")
                ) <= active_tolerance:
                    raise ReceiptError(f"{study}.{label}[{index}] stores an inactive nonzero flux")
            if set(active) != nonzero_ids:
                raise ReceiptError(f"{study}.{label}[{index}].active_by_flux disagrees with fluxes")
            union.update(active)
        declared = set(_strings(corneto.get(f"{label}_active_union"), f"{study}.{label}_active_union"))
        if declared != union:
            raise ReceiptError(f"{study}: declared {label} active union disagrees with summaries")
        if _integer(corneto.get(f"{label}_active_union_size"), f"{study}.{label}_active_union_size") != len(union):
            raise ReceiptError(f"{study}: declared {label} active union size disagrees")
        active_unions[label] = union
        status_counts[label] = dict(sorted(counts.items()))

    return {
        "study": study,
        "path": str(path),
        "sample_count": len(sample_ids),
        "run_accessions": sample_ids,
        "solver_requested": requested,
        "solver_used": used,
        "candidate_count": len(candidate_ids),
        "candidate_ids": sorted(candidate_ids),
        "status_counts": status_counts,
        "independent_active_union": sorted(active_unions["independent"]),
        "joint_active_union": sorted(active_unions["joint"]),
    }


def validate(specs: list[str]) -> dict[str, Any]:
    if not specs:
        raise ReceiptError("at least one --receipt is required")
    parsed = [_parse_spec(spec) for spec in specs]
    studies = [study for study, _ in parsed]
    if len(studies) != len(set(studies)):
        raise ReceiptError("duplicate study in --receipt")
    if set(studies) != set(EXPECTED_SAMPLE_COUNTS):
        raise ReceiptError(
            f"studies must exactly match the frozen four-cohort contract: {sorted(EXPECTED_SAMPLE_COUNTS)}"
        )
    summaries = [_read_receipt(study, path) for study, path in parsed]
    pairwise: dict[str, dict[str, float | int | None]] = {}
    for index, left in enumerate(summaries):
        for right in summaries[index + 1 :]:
            key = f"{left['study']}__{right['study']}"
            left_ids = set(left["candidate_ids"])
            right_ids = set(right["candidate_ids"])
            left_flux = set(left["joint_active_union"])
            right_flux = set(right["joint_active_union"])
            pairwise[key] = {
                "candidate_intersection": len(left_ids & right_ids),
                "candidate_union": len(left_ids | right_ids),
                "candidate_jaccard": _jaccard(left_ids, right_ids),
                "joint_active_intersection": len(left_flux & right_flux),
                "joint_active_union": len(left_flux | right_flux),
                "joint_active_jaccard": _jaccard(left_flux, right_flux),
            }
    return {
        "status": "valid",
        "response_blind": True,
        "cohort_count": len(summaries),
        "total_sample_count": sum(item["sample_count"] for item in summaries),
        "contract": {
            "expected_sample_counts": EXPECTED_SAMPLE_COUNTS,
            "expected_solver": EXPECTED_SOLVER,
            "max_candidates": EXPECTED_CANDIDATES,
            "growth_fraction": EXPECTED_GROWTH_FRACTION,
            "independent_lambda": EXPECTED_INDEPENDENT_LAMBDA,
            "joint_lambda": EXPECTED_JOINT_LAMBDA,
        },
        "cohorts": summaries,
        "pairwise": pairwise,
        "claim_limit": "descriptive metabolic state and solver receipt comparison; no Taxol phenotype or causal claim",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", action="append", required=True, metavar="STUDY=PATH")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = validate(args.receipt)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.output.with_suffix(args.output.suffix + ".tmp")
        tmp.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        tmp.replace(args.output)
    except ReceiptError as error:
        print(json.dumps({"status": "invalid", "error": str(error)}, indent=2), file=sys.stderr)
        return 1
    print(json.dumps({"status": result["status"], "cohort_count": result["cohort_count"], "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
