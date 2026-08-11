#!/usr/bin/env python3
"""Validate and compare pooled and cohort-joint CARNIVAL lambda grids."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

EXPECTED_STUDIES = (
    "E-MTAB-7223",
    "E-MTAB-10801",
    "E-MTAB-11000",
    "E-MTAB-14568",
)
EXPECTED_LAMBDAS = (0.0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0)
Edge = tuple[str, str, int]


def _receipt_spec(value: str) -> tuple[str, str | None, float, Path]:
    pieces = value.split("|", 3)
    if len(pieces) != 4:
        raise argparse.ArgumentTypeError("receipt must be MODE|STUDY|LAMBDA|PATH")
    mode, raw_study, raw_lambda, raw_path = pieces
    if mode not in {"pooled", "cohort"}:
        raise argparse.ArgumentTypeError("MODE must be pooled or cohort")
    study = raw_study or None
    if mode == "pooled" and study is not None:
        raise argparse.ArgumentTypeError("pooled receipt STUDY must be empty")
    if mode == "cohort" and study not in EXPECTED_STUDIES:
        raise argparse.ArgumentTypeError(f"unknown cohort study: {study}")
    return mode, study, float(raw_lambda), Path(raw_path)


def _edges(condition: dict[str, Any]) -> set[Edge]:
    return {
        (str(edge["source"]), str(edge["target"]), int(edge.get("sign", 0)))
        for edge in condition.get("selected_edges", [])
    }


def _jaccard(left: set[Edge], right: set[Edge]) -> float:
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def _condition_map(receipt: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = receipt.get("conditions")
    if not isinstance(rows, list) or not rows:
        raise ValueError("receipt has no conditions")
    result = {str(row["run_accession"]): row for row in rows}
    if len(result) != len(rows):
        raise ValueError("receipt has duplicate conditions")
    return result


def summarize(
    specs: list[tuple[str, str | None, float, Path]], output: Path
) -> dict[str, Any]:
    if output.exists():
        raise ValueError(f"refusing to overwrite {output}")
    receipts: dict[tuple[str, str | None, float], dict[str, Any]] = {}
    bundle_hashes: set[str] = set()
    for mode, study, lambda_value, path in specs:
        key = (mode, study, lambda_value)
        if key in receipts:
            raise ValueError(f"duplicate receipt: {key}")
        data = json.loads(path.read_text(encoding="utf-8"))
        if data.get("schema_version") != "regulatory_multisample_solution.v1":
            raise ValueError(f"{path}: unexpected schema")
        if data.get("status") != "completed":
            raise ValueError(f"{path}: status={data.get('status')!r}")
        method = data.get("method", {})
        if not method.get("single_joint_problem"):
            raise ValueError(f"{path}: not a joint problem")
        if method.get("lambda_scaling") != "mean_fit":
            raise ValueError(f"{path}: lambda scaling is not mean_fit")
        if abs(float(method.get("lambda_nominal")) - lambda_value) > 1e-12:
            raise ValueError(f"{path}: nominal lambda mismatch")
        expected_solver_lambda = lambda_value * int(method.get("condition_count"))
        if abs(float(method.get("lambda_solver")) - expected_solver_lambda) > 1e-9:
            raise ValueError(f"{path}: solver lambda is not condition-normalized")
        bundle_hashes.add(str(data.get("bundle", {}).get("sha256")))
        data["_path"] = str(path)
        data["_condition_map"] = _condition_map(data)
        receipts[key] = data
    if len(bundle_hashes) != 1 or "None" in bundle_hashes:
        raise ValueError(f"receipts do not share one frozen bundle: {sorted(bundle_hashes)}")

    expected_keys = {
        (mode, study, lambda_value)
        for mode, studies in (("pooled", (None,)), ("cohort", EXPECTED_STUDIES))
        for study in studies
        for lambda_value in EXPECTED_LAMBDAS
    }
    missing = expected_keys - set(receipts)
    extra = set(receipts) - expected_keys
    if missing or extra:
        raise ValueError(f"grid mismatch; missing={sorted(missing)}, extra={sorted(extra)}")

    pooled_runs_by_lambda = {
        value: set(receipts[("pooled", None, value)]["_condition_map"])
        for value in EXPECTED_LAMBDAS
    }
    if len({frozenset(runs) for runs in pooled_runs_by_lambda.values()}) != 1:
        raise ValueError("pooled condition set changes across lambda")
    for study in EXPECTED_STUDIES:
        run_sets = {
            value: set(receipts[("cohort", study, value)]["_condition_map"])
            for value in EXPECTED_LAMBDAS
        }
        if len({frozenset(runs) for runs in run_sets.values()}) != 1:
            raise ValueError(f"{study}: cohort condition set changes across lambda")
        if set().union(*run_sets.values()) - pooled_runs_by_lambda[0.0]:
            raise ValueError(f"{study}: cohort runs absent from pooled receipt")

    modes: dict[str, Any] = {}
    for mode, study in (("pooled", None), *( ("cohort", value) for value in EXPECTED_STUDIES)):
        label = "pooled" if mode == "pooled" else str(study)
        baseline = receipts[(mode, study, 0.0)]
        baseline_conditions = baseline["_condition_map"]
        baseline_union = set().union(*(_edges(row) for row in baseline_conditions.values()))
        rows: dict[str, Any] = {}
        for lambda_value in EXPECTED_LAMBDAS:
            receipt = receipts[(mode, study, lambda_value)]
            conditions = receipt["_condition_map"]
            union = set().union(*(_edges(row) for row in conditions.values()))
            per_sample = [
                _jaccard(_edges(baseline_conditions[run]), _edges(conditions[run]))
                for run in sorted(conditions)
            ]
            agreements = [
                float(row["output_sign_agreement"])
                for row in conditions.values()
                if row.get("output_sign_agreement") is not None
            ]
            rows[str(lambda_value)] = {
                "path": receipt["_path"],
                "solver_status": receipt.get("solver", {}).get("status"),
                "included_conditions": receipt.get("scope_counts", {}).get(
                    "included_conditions"
                ),
                "lambda_solver": receipt.get("method", {}).get("lambda_solver"),
                "edge_union_size": len(union),
                "edge_union_jaccard_vs_lambda0": _jaccard(baseline_union, union),
                "mean_sample_edge_jaccard_vs_lambda0": sum(per_sample) / len(per_sample),
                "mean_output_sign_agreement": (
                    sum(agreements) / len(agreements) if agreements else None
                ),
                "condition_status_counts": dict(
                    Counter(str(row.get("status")) for row in conditions.values())
                ),
            }
        modes[label] = rows

    pooled_vs_cohort: dict[str, Any] = {}
    for lambda_value in EXPECTED_LAMBDAS:
        pooled = receipts[("pooled", None, lambda_value)]["_condition_map"]
        by_study: dict[str, Any] = {}
        all_jaccards: list[float] = []
        for study in EXPECTED_STUDIES:
            cohort = receipts[("cohort", study, lambda_value)]["_condition_map"]
            values = [
                _jaccard(_edges(pooled[run]), _edges(cohort[run]))
                for run in sorted(cohort)
            ]
            all_jaccards.extend(values)
            by_study[study] = {
                "condition_count": len(values),
                "mean_sample_edge_jaccard": sum(values) / len(values),
                "sample_edge_jaccards": dict(zip(sorted(cohort), values, strict=True)),
            }
        pooled_vs_cohort[str(lambda_value)] = {
            "mean_sample_edge_jaccard": sum(all_jaccards) / len(all_jaccards),
            "by_study": by_study,
        }

    result = {
        "schema_version": "regulatory_multisample_grid_summary.v1",
        "status": "completed",
        "response_blind": True,
        "bundle_sha256": next(iter(bundle_hashes)),
        "lambda_values": list(EXPECTED_LAMBDAS),
        "lambda_scaling": "mean condition fit + nominal lambda * union edge count",
        "modes": modes,
        "pooled_vs_cohort": pooled_vs_cohort,
        "claim_limit": (
            "method stability only; lambda must not be selected using held-out phenotype"
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--receipt", action="append", type=_receipt_spec, required=True
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.receipt, args.output)
    print(
        json.dumps(
            {
                "status": result["status"],
                "lambda_count": len(result["lambda_values"]),
                "mode_count": len(result["modes"]),
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
