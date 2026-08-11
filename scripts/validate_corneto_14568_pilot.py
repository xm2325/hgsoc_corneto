#!/usr/bin/env python3
"""Validate a CORNETO 14568 pilot receipt without rerunning a solver.

The validator is intentionally standard-library-only and reads one JSON
receipt.  It checks the invariants recorded by ``run_corneto_14568_pilot.py``
but never opens a model, reads a solver license, or prints sample rows.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


class ReceiptValidationError(Exception):
    """Raised for malformed receipt structures during validation."""


def _mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ReceiptValidationError(f"{name} must be an object")
    return value


def _sequence(value: Any, name: str) -> list[Any]:
    if not isinstance(value, list):
        raise ReceiptValidationError(f"{name} must be an array")
    return value


def _number(value: Any, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool):
        raise ReceiptValidationError(f"{name} must be numeric")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ReceiptValidationError(f"{name} must be numeric") from error
    if not math.isfinite(parsed):
        raise ReceiptValidationError(f"{name} must be finite")
    if positive and parsed <= 0:
        raise ReceiptValidationError(f"{name} must be positive")
    return parsed


def _integer(value: Any, name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ReceiptValidationError(f"{name} must be an integer")
    if value < minimum:
        raise ReceiptValidationError(f"{name} must be >= {minimum}")
    return value


def _string_list(value: Any, name: str, *, require_unique: bool = True) -> list[str]:
    values = _sequence(value, name)
    if any(not isinstance(item, str) or not item for item in values):
        raise ReceiptValidationError(f"{name} must contain non-empty strings")
    result = list(values)
    if require_unique and len(result) != len(set(result)):
        raise ReceiptValidationError(f"{name} contains duplicate IDs")
    return result


def _same_keys(actual: dict[str, Any], expected: set[str], name: str, errors: list[str]) -> None:
    keys = set(actual)
    missing = sorted(expected - keys)
    extra = sorted(keys - expected)
    if missing:
        errors.append(f"{name} missing keys: {missing}")
    if extra:
        errors.append(f"{name} has unexpected keys: {extra}")


def _check_summary(
    summary: Any,
    *,
    name: str,
    expected_condition: str,
    active_tolerance: float,
    errors: list[str],
) -> set[str]:
    try:
        item = _mapping(summary, name)
        condition = item.get("condition")
        if condition != expected_condition:
            errors.append(f"{name}.condition={condition!r}, expected {expected_condition!r}")
        status = str(item.get("status", "")).casefold()
        if status not in {"optimal", "optimal_inaccurate"}:
            errors.append(f"{name}.status={status!r}")
        _number(item.get("problem_objective_value"), f"{name}.problem_objective_value")

        active_by_flux = set(
            _string_list(item.get("active_by_flux"), f"{name}.active_by_flux")
        )
        _string_list(item.get("active_by_indicator"), f"{name}.active_by_indicator")
        nonzero = _sequence(item.get("nonzero_fluxes"), f"{name}.nonzero_fluxes")
        nonzero_ids: set[str] = set()
        for index, pair in enumerate(nonzero):
            if not isinstance(pair, list) or len(pair) != 2:
                raise ReceiptValidationError(f"{name}.nonzero_fluxes[{index}] must be [id, value]")
            reaction_id, flux = pair
            if not isinstance(reaction_id, str) or not reaction_id:
                raise ReceiptValidationError(f"{name}.nonzero_fluxes[{index}] has invalid ID")
            if reaction_id in nonzero_ids:
                raise ReceiptValidationError(f"{name}.nonzero_fluxes contains duplicate IDs")
            nonzero_ids.add(reaction_id)
            flux_value = _number(flux, f"{name}.nonzero_fluxes[{index}][1]")
            if abs(flux_value) <= active_tolerance:
                raise ReceiptValidationError(
                    f"{name}.nonzero_fluxes[{index}] is not above active tolerance"
                )
        if nonzero_ids != active_by_flux:
            errors.append(f"{name}.active_by_flux disagrees with nonzero_fluxes IDs")
        return active_by_flux
    except ReceiptValidationError as error:
        errors.append(str(error))
        return set()


def validate(
    receipt: dict[str, Any],
    *,
    expected_solver: str,
    expected_samples: int,
    expected_candidates: int,
    expected_growth_fraction: float,
) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    if receipt.get("status") != "completed":
        errors.append(f"status={receipt.get('status')!r}, expected 'completed'")
    if receipt.get("study_accession") != "E-MTAB-14568":
        errors.append(f"study_accession={receipt.get('study_accession')!r}")
    if receipt.get("primary_only") is not True:
        errors.append("primary_only must be true")

    try:
        sample_count = _integer(receipt.get("sample_count"), "sample_count", minimum=1)
    except ReceiptValidationError as error:
        errors.append(str(error))
        sample_count = 0
    if sample_count != expected_samples:
        errors.append(f"sample_count={sample_count}, expected {expected_samples}")

    try:
        samples = _sequence(receipt.get("samples"), "samples")
        if len(samples) != sample_count:
            errors.append("samples length disagrees with sample_count")
        run_ids: list[str] = []
        for index, row in enumerate(samples):
            item = _mapping(row, f"samples[{index}]")
            run_id = item.get("run_accession")
            if not isinstance(run_id, str) or not run_id:
                raise ReceiptValidationError(f"samples[{index}].run_accession is invalid")
            run_ids.append(run_id)
        if len(run_ids) != len(set(run_ids)):
            errors.append("samples contain duplicate run_accession values")
    except ReceiptValidationError as error:
        errors.append(str(error))
        run_ids = []

    try:
        solver = _mapping(receipt.get("solver"), "solver")
        requested = str(solver.get("requested", "")).casefold()
        used = str(solver.get("used", "")).casefold()
        expected_solver = expected_solver.casefold()
        if requested != expected_solver:
            errors.append(f"solver.requested={requested!r}, expected {expected_solver!r}")
        if used != expected_solver:
            errors.append(f"solver.used={used!r}, expected {expected_solver!r}")
        available = [str(value).casefold() for value in _sequence(solver.get("available"), "solver.available")]
        if expected_solver not in available:
            errors.append(f"solver.available does not contain {expected_solver!r}")
        if solver.get("fallback_reason") is not None:
            errors.append("solver.fallback_reason is not null")
        if solver.get("solve_fallback") is not None:
            errors.append("solver.solve_fallback is not null")
    except ReceiptValidationError as error:
        errors.append(str(error))

    try:
        candidate = _mapping(receipt.get("candidate_selection"), "candidate_selection")
        max_candidates = _integer(candidate.get("max_candidates"), "candidate_selection.max_candidates")
        if max_candidates != expected_candidates:
            errors.append(
                f"candidate_selection.max_candidates={max_candidates}, expected {expected_candidates}"
            )
        selected = _string_list(
            candidate.get("selected_reaction_ids"),
            "candidate_selection.selected_reaction_ids",
        )
        selected_count = _integer(candidate.get("selected_count"), "candidate_selection.selected_count", minimum=1)
        if selected_count != len(selected):
            errors.append("selected_count disagrees with selected_reaction_ids")
        if selected_count > expected_candidates:
            errors.append("selected_count exceeds max candidate budget")
        selected_scores = _mapping(
            candidate.get("selected_median_proposed_upper"),
            "candidate_selection.selected_median_proposed_upper",
        )
        _same_keys(selected_scores, set(selected), "selected_median_proposed_upper", errors)
        for reaction_id in selected:
            try:
                _number(selected_scores.get(reaction_id), f"selected_median_proposed_upper[{reaction_id}]", positive=True)
            except ReceiptValidationError as error:
                errors.append(str(error))
        audits = _mapping(candidate.get("audits_by_sample"), "candidate_selection.audits_by_sample")
        _same_keys(audits, set(run_ids), "audits_by_sample", errors)
        for run_id in run_ids:
            try:
                audit = _mapping(audits.get(run_id), f"audits_by_sample[{run_id}]")
                all_count = _integer(audit.get("all_candidates"), f"audits_by_sample[{run_id}].all_candidates")
                positive_count = _integer(audit.get("positive_candidates"), f"audits_by_sample[{run_id}].positive_candidates")
                if positive_count > all_count:
                    errors.append(f"audits_by_sample[{run_id}] positive count exceeds all count")
            except ReceiptValidationError as error:
                errors.append(str(error))
        if candidate.get("bounds_clamped_to_human_gem") is not True:
            errors.append("bounds_clamped_to_human_gem must be true")
        skipped = _mapping(
            candidate.get("bounds_skipped_as_empty_intersection"),
            "bounds_skipped_as_empty_intersection",
        )
        _same_keys(skipped, set(run_ids), "bounds_skipped_as_empty_intersection", errors)
        for run_id in run_ids:
            try:
                _integer(skipped.get(run_id), f"bounds_skipped_as_empty_intersection[{run_id}]")
            except ReceiptValidationError as error:
                errors.append(str(error))
        warnings.append(
            "Human-GEM original bounds are not stored in this receipt; the clamp check verifies the runner's declared invariant and skip accounting."
        )
    except ReceiptValidationError as error:
        errors.append(str(error))

    try:
        objective = _mapping(receipt.get("objective"), "objective")
        coefficient = _number(objective.get("biomass_coefficient"), "objective.biomass_coefficient")
        if coefficient != -1.0:
            errors.append(f"objective.biomass_coefficient={coefficient}, expected -1.0")
        growth_fraction = _number(objective.get("growth_fraction"), "objective.growth_fraction")
        if not math.isclose(growth_fraction, expected_growth_fraction, rel_tol=0, abs_tol=1e-12):
            errors.append(f"objective.growth_fraction={growth_fraction}, expected {expected_growth_fraction}")
        optima = _mapping(objective.get("growth_optima"), "objective.growth_optima")
        targets = _mapping(objective.get("growth_targets"), "objective.growth_targets")
        _same_keys(optima, set(run_ids), "growth_optima", errors)
        _same_keys(targets, set(run_ids), "growth_targets", errors)
        for run_id in run_ids:
            try:
                optimum = _number(optima.get(run_id), f"growth_optima[{run_id}]", positive=True)
                target = _number(targets.get(run_id), f"growth_targets[{run_id}]", positive=True)
                if not math.isclose(target, expected_growth_fraction * optimum, rel_tol=1e-9, abs_tol=1e-10):
                    errors.append(f"growth_targets[{run_id}] is not growth_fraction * growth_optima")
            except ReceiptValidationError as error:
                errors.append(str(error))
        for field in ("independent_lambda", "joint_lambda"):
            try:
                _number(objective.get(field), f"objective.{field}")
                if float(objective[field]) < 0:
                    errors.append(f"objective.{field} must be non-negative")
            except ReceiptValidationError as error:
                errors.append(str(error))
    except ReceiptValidationError as error:
        errors.append(str(error))

    active_unions: dict[str, set[str]] = {"independent": set(), "joint": set()}
    try:
        corneto = _mapping(receipt.get("corneto"), "corneto")
        if str(corneto.get("solver", "")).casefold() != expected_solver:
            errors.append("corneto.solver disagrees with expected solver")
        conditions = _string_list(corneto.get("conditions"), "corneto.conditions")
        if conditions != run_ids:
            errors.append("corneto.conditions disagrees with sample run order")
        active_tolerance = _number(corneto.get("active_tolerance"), "corneto.active_tolerance", positive=True)
        for label in ("independent", "joint"):
            summaries = _sequence(corneto.get(label), f"corneto.{label}")
            expected_length = len(run_ids)
            if len(summaries) != expected_length:
                errors.append(f"corneto.{label} length={len(summaries)}, expected {expected_length}")
            for index, run_id in enumerate(run_ids[: len(summaries)]):
                active_unions[label].update(
                    _check_summary(
                        summaries[index],
                        name=f"corneto.{label}[{index}]",
                        expected_condition=run_id,
                        active_tolerance=active_tolerance,
                        errors=errors,
                    )
                )
        for label in ("independent", "joint"):
            declared = _string_list(corneto.get(f"{label}_active_union"), f"corneto.{label}_active_union")
            if set(declared) != active_unions[label] or declared != sorted(declared):
                errors.append(f"corneto.{label}_active_union disagrees with condition summaries")
            declared_size = _integer(
                corneto.get(f"{label}_active_union_size"),
                f"corneto.{label}_active_union_size",
            )
            if declared_size != len(declared):
                errors.append(f"corneto.{label}_active_union_size disagrees with union")
    except ReceiptValidationError as error:
        errors.append(str(error))

    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "checks": {
            "expected_solver": expected_solver,
            "expected_samples": expected_samples,
            "expected_candidates": expected_candidates,
            "expected_growth_fraction": expected_growth_fraction,
            "run_accessions": run_ids,
            "independent_active_union_size": len(active_unions["independent"]),
            "joint_active_union_size": len(active_unions["joint"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--expected-solver", default="highs")
    parser.add_argument("--expected-samples", type=int, default=1)
    parser.add_argument("--expected-candidates", type=int, default=25)
    parser.add_argument("--expected-growth-fraction", type=float, default=0.9)
    args = parser.parse_args()
    try:
        receipt = json.loads(args.receipt.read_text(encoding="utf-8"))
        receipt_mapping = _mapping(receipt, "receipt")
        result = validate(
            receipt_mapping,
            expected_solver=args.expected_solver,
            expected_samples=args.expected_samples,
            expected_candidates=args.expected_candidates,
            expected_growth_fraction=args.expected_growth_fraction,
        )
    except (OSError, json.JSONDecodeError, ReceiptValidationError) as error:
        result = {"valid": False, "errors": [str(error)], "warnings": [], "checks": {}}
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    sys.exit(main())

