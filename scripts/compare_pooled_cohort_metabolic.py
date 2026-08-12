#!/usr/bin/env python3
"""Compare pooled-60 joint metabolic CORNETO with four cohort-joint receipts."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


EXPECTED = {
    "E-MTAB-7223": 9,
    "E-MTAB-10801": 13,
    "E-MTAB-11000": 11,
    "E-MTAB-14568": 27,
}


def _object(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _strings(value: Any, label: str) -> list[str]:
    if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
        raise ValueError(f"{label} must be an array of non-empty strings")
    if len(value) != len(set(value)):
        raise ValueError(f"{label} contains duplicates")
    return list(value)


def _jaccard(left: set[str], right: set[str]) -> float | None:
    union = left | right
    return None if not union else len(left & right) / len(union)


def _parse_spec(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"cohort receipt must be STUDY=PATH, got {value!r}")
    study, raw_path = value.split("=", 1)
    if study not in EXPECTED or not raw_path:
        raise ValueError(f"unsupported cohort receipt {value!r}")
    return study, Path(raw_path)


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        return _object(json.loads(path.read_text(encoding="utf-8")), label)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read {label}: {error}") from error


def _joint_by_condition(root: dict[str, Any], label: str) -> dict[str, set[str]]:
    corneto = _object(root.get("corneto"), f"{label}.corneto")
    summaries = corneto.get("joint")
    if not isinstance(summaries, list):
        raise ValueError(f"{label}.corneto.joint must be an array")
    result: dict[str, set[str]] = {}
    for index, raw in enumerate(summaries):
        item = _object(raw, f"{label}.joint[{index}]")
        condition = item.get("condition")
        status = str(item.get("status", "")).casefold()
        if not isinstance(condition, str) or not condition or condition in result:
            raise ValueError(f"{label}.joint[{index}] has an invalid condition")
        if status not in {"optimal", "optimal_inaccurate"}:
            raise ValueError(f"{label}.{condition} is not optimal: {status!r}")
        result[condition] = set(_strings(item.get("active_by_flux"), f"{label}.{condition}.active"))
    return result


def compare(pooled_path: Path, specs: list[str]) -> dict[str, Any]:
    parsed = [_parse_spec(value) for value in specs]
    if len(parsed) != 4 or {study for study, _ in parsed} != set(EXPECTED):
        raise ValueError("cohort receipts must name all four frozen studies exactly once")
    pooled = _read_json(pooled_path, "pooled receipt")
    if pooled.get("status") != "completed" or pooled.get("joint_only") is not True:
        raise ValueError("pooled receipt is not a completed joint-only analysis")
    if pooled.get("sample_count") != 60:
        raise ValueError("pooled receipt does not contain exactly 60 conditions")
    pooled_solver = _object(pooled.get("solver"), "pooled.solver")
    if str(pooled_solver.get("used", "")).casefold() != "gurobi":
        raise ValueError("pooled receipt did not use Gurobi")
    pooled_candidate = set(
        _strings(
            _object(pooled.get("candidate_selection"), "pooled.candidates").get(
                "selected_reaction_ids"
            ),
            "pooled selected reactions",
        )
    )
    pooled_joint = _joint_by_condition(pooled, "pooled")
    if len(pooled_joint) != 60:
        raise ValueError("pooled joint summaries do not contain 60 unique conditions")

    cohorts: dict[str, Any] = {}
    all_cohort_runs: set[str] = set()
    for study, path in parsed:
        root = _read_json(path, f"{study} receipt")
        if root.get("status") != "completed" or root.get("study_accession") != study:
            raise ValueError(f"{study} receipt is not completed or is mislabeled")
        if root.get("sample_count") != EXPECTED[study]:
            raise ValueError(f"{study} sample count differs from the frozen contract")
        candidate = set(
            _strings(
                _object(root.get("candidate_selection"), f"{study}.candidates").get(
                    "selected_reaction_ids"
                ),
                f"{study} selected reactions",
            )
        )
        joint = _joint_by_condition(root, study)
        if len(joint) != EXPECTED[study]:
            raise ValueError(f"{study} joint condition count is inconsistent")
        overlap = set(joint) & all_cohort_runs
        if overlap:
            raise ValueError(f"run accessions occur in multiple cohort receipts: {sorted(overlap)}")
        all_cohort_runs.update(joint)
        per_sample = {
            run_id: _jaccard(active, pooled_joint[run_id])
            for run_id, active in joint.items()
            if run_id in pooled_joint
        }
        if len(per_sample) != EXPECTED[study]:
            raise ValueError(f"{study} runs do not all occur in the pooled receipt")
        finite = [value for value in per_sample.values() if value is not None and math.isfinite(value)]
        cohort_union = set().union(*joint.values()) if joint else set()
        pooled_study_union = set().union(*(pooled_joint[run] for run in joint)) if joint else set()
        cohorts[study] = {
            "sample_count": len(joint),
            "candidate_count": len(candidate),
            "pooled_candidate_count": len(pooled_candidate),
            "candidate_jaccard": _jaccard(candidate, pooled_candidate),
            "cohort_joint_union_size": len(cohort_union),
            "pooled_same_samples_union_size": len(pooled_study_union),
            "joint_union_jaccard": _jaccard(cohort_union, pooled_study_union),
            "mean_sample_active_jaccard": sum(finite) / len(finite) if finite else None,
            "per_sample_active_jaccard": per_sample,
        }
    if all_cohort_runs != set(pooled_joint):
        raise ValueError("four cohort run union does not exactly match pooled-60 conditions")
    return {
        "status": "completed",
        "response_blind": True,
        "comparison": "pooled_60_joint_vs_four_cohort_joint_metabolic_sparse_fba",
        "sample_count": 60,
        "cohorts": cohorts,
        "claim_limit": (
            "Internal robustness comparison of models derived from the same RNA inputs. "
            "Agreement is not independent validation, measured flux, or phenotype evidence."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pooled", type=Path, required=True)
    parser.add_argument("--cohort", action="append", required=True, metavar="STUDY=PATH")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        print(f"refusing to overwrite output: {args.output}", file=sys.stderr)
        return 1
    try:
        result = compare(args.pooled, args.cohort)
    except ValueError as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"status": "completed", "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
