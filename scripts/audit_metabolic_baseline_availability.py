#!/usr/bin/env python3
"""Write a fail-transparent snapshot of four metabolic baseline receipts.

Unlike the strict final comparator, this audit always writes a receipt.  A
missing upstream result is recorded as ``incomplete`` rather than silently
turning a scheduler dependency into a scientific success.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


EXPECTED = {
    "E-MTAB-7223": 9,
    "E-MTAB-10801": 13,
    "E-MTAB-11000": 11,
    "E-MTAB-14568": 27,
}


def _parse_spec(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"receipt must be STUDY=PATH, got {value!r}")
    study, raw_path = value.split("=", 1)
    if study not in EXPECTED or not raw_path:
        raise ValueError(f"unsupported or incomplete receipt spec {value!r}")
    return study, Path(raw_path)


def _inspect(study: str, path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "study_accession": study,
        "path": str(path),
        "expected_sample_count": EXPECTED[study],
    }
    if not path.is_file():
        result.update(state="missing", valid_for_final_comparison=False)
        return result
    try:
        root = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        result.update(
            state="invalid_json",
            valid_for_final_comparison=False,
            error=f"{type(error).__name__}: {error}",
        )
        return result
    if not isinstance(root, dict):
        result.update(state="invalid_schema", valid_for_final_comparison=False)
        return result
    solver = root.get("solver") if isinstance(root.get("solver"), dict) else {}
    candidate = (
        root.get("candidate_selection")
        if isinstance(root.get("candidate_selection"), dict)
        else {}
    )
    objective = root.get("objective") if isinstance(root.get("objective"), dict) else {}
    checks = {
        "status_completed": root.get("status") == "completed",
        "study_matches": root.get("study_accession") == study,
        "sample_count_matches": root.get("sample_count") == EXPECTED[study],
        "primary_only": root.get("primary_only") is True,
        "solver_requested_gurobi": str(solver.get("requested", "")).casefold()
        == "gurobi",
        "solver_used_gurobi": str(solver.get("used", "")).casefold() == "gurobi",
        "no_solver_fallback": solver.get("fallback_reason") is None
        and solver.get("solve_fallback") is None,
        "candidate_budget_25": candidate.get("max_candidates") == 25,
        "growth_fraction_0p9": objective.get("growth_fraction") == 0.9,
        "independent_lambda_0p1": objective.get("independent_lambda") == 0.1,
        "joint_lambda_1": objective.get("joint_lambda") == 1.0,
    }
    valid = all(checks.values())
    result.update(
        state="ready" if valid else "invalid_contract",
        valid_for_final_comparison=valid,
        checks=checks,
        observed={
            "status": root.get("status"),
            "sample_count": root.get("sample_count"),
            "solver_requested": solver.get("requested"),
            "solver_used": solver.get("used"),
            "selected_count": candidate.get("selected_count"),
        },
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", action="append", required=True, metavar="STUDY=PATH")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        specs = [_parse_spec(value) for value in args.receipt]
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 1
    studies = [study for study, _ in specs]
    if len(studies) != len(set(studies)) or set(studies) != set(EXPECTED):
        print("receipts must name each frozen cohort exactly once", file=sys.stderr)
        return 1
    cohorts = [_inspect(study, path) for study, path in specs]
    ready = sum(bool(row["valid_for_final_comparison"]) for row in cohorts)
    invalid = sum(row["state"].startswith("invalid") for row in cohorts)
    status = "complete" if ready == len(EXPECTED) else ("invalid" if invalid else "incomplete")
    result = {
        "status": status,
        "response_blind": True,
        "ready_cohort_count": ready,
        "expected_cohort_count": len(EXPECTED),
        "cohorts": cohorts,
        "final_comparison_permitted": ready == len(EXPECTED),
        "claim_limit": (
            "Scheduler-independent availability audit only. Biological cross-cohort "
            "claims require all four strict receipts to validate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"status": status, "ready": ready, "output": str(args.output)}))
    return 0 if status == "complete" else (1 if status == "invalid" else 2)


if __name__ == "__main__":
    sys.exit(main())
