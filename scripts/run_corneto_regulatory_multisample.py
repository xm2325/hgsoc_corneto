#!/usr/bin/env python3
"""Solve one true multi-condition CORNETO CarnivalFlow problem.

Unlike the legacy independent pilot, this script builds one ``Data`` object
containing all selected conditions and calls ``CarnivalFlow`` once.  The
structured regularizer therefore acts on the union of edge signals across
conditions, as intended by multi-sample CARNIVAL.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if hasattr(value, "tolist"):
        return value.tolist()
    try:
        return float(value)
    except (TypeError, ValueError):
        return str(value)


def _effective_lambda(nominal_lambda: float, condition_count: int, scaling: str) -> float:
    if nominal_lambda < 0:
        raise ValueError("lambda must be non-negative")
    if condition_count < 1:
        raise ValueError("condition_count must be positive")
    if scaling == "mean_fit":
        return nominal_lambda * condition_count
    if scaling == "raw":
        return nominal_lambda
    raise ValueError(f"unknown lambda scaling: {scaling}")


def _edge_records(
    method: Any, edge_matrix: Any, condition_names: list[str]
) -> dict[str, list[dict[str, object]]]:
    import numpy as np

    values = np.asarray(edge_matrix, dtype=float)
    edge_count = len(method.processed_graph.E)
    if values.ndim == 1:
        values = values.reshape((-1, 1))
    if values.shape == (len(condition_names), edge_count):
        values = values.T
    if values.shape != (edge_count, len(condition_names)):
        raise ValueError(
            f"edge_has_signal shape {values.shape} != ({edge_count}, {len(condition_names)})"
        )
    selected: dict[str, list[dict[str, object]]] = {name: [] for name in condition_names}
    for edge_index in range(edge_count):
        source_set, target_set = method.processed_graph.E[edge_index]
        sources, targets = sorted(source_set), sorted(target_set)
        if not sources or not targets:
            continue
        sign = method.processed_graph.get_attr_edge(edge_index).get("interaction")
        if sign is not None:
            sign = int(sign)
        for condition_index, name in enumerate(condition_names):
            if values[edge_index, condition_index] > 0.5:
                selected[name].append(
                    {"source": sources[0], "target": targets[0], "sign": sign}
                )
    return selected


def solve_bundle(
    *,
    bundle_path: Path,
    output: Path,
    mode: str,
    study: str | None,
    nominal_lambda: float,
    lambda_scaling: str,
    solver: str,
    max_flow: int,
    max_seconds: int,
    max_conditions: int | None,
) -> dict[str, object]:
    if output.exists():
        raise ValueError(f"refusing to overwrite {output}")
    bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
    if bundle.get("schema_version") != "regulatory_multisample_input.v1":
        raise ValueError("unsupported or missing multi-sample input schema")
    all_conditions = bundle.get("conditions")
    graph_rows = bundle.get("graph")
    if not isinstance(all_conditions, list) or not isinstance(graph_rows, list):
        raise ValueError("bundle is missing conditions or graph")
    if mode == "cohort" and not study:
        raise ValueError("--study is required for cohort mode")
    if mode == "pooled" and study:
        raise ValueError("--study is not allowed for pooled mode")
    scoped = [
        row
        for row in all_conditions
        if mode == "pooled" or row.get("study_accession") == study
    ]
    if max_conditions is not None:
        scoped = scoped[:max_conditions]
    included = [row for row in scoped if row.get("preprocessing_status") == "included"]
    if len(included) < 2:
        raise ValueError(
            "true multi-sample solve requires at least two included conditions; "
            f"found {len(included)}"
        )
    condition_names = [str(row["run_accession"]) for row in included]
    if len(condition_names) != len(set(condition_names)):
        raise ValueError("duplicate run_accession in selected conditions")

    # CORNETO 1.0 currently adds one unweighted error objective per condition
    # and lambda_solver * |edge union|.  To implement the paper-scale objective
    # mean(error) + lambda_nominal * |edge union|, multiply lambda by n.
    solver_lambda = _effective_lambda(nominal_lambda, len(included), lambda_scaling)
    if lambda_scaling == "mean_fit":
        scaling_formula = "lambda_solver = lambda_nominal * included_condition_count"
    else:
        solver_lambda = nominal_lambda
        scaling_formula = "lambda_solver = lambda_nominal (raw CORNETO objective scale)"

    started = datetime.now(UTC).replace(microsecond=0).isoformat()
    receipt: dict[str, object] = {
        "schema_version": "regulatory_multisample_solution.v1",
        "status": "running",
        "started_at_utc": started,
        "response_blind": True,
        "analysis_mode": mode,
        "study_accession": study,
        "bundle": {"path": str(bundle_path), "sha256": _sha(bundle_path)},
        "scope_counts": {
            "manifest_conditions": len(scoped),
            "included_conditions": len(included),
            "preprocessing_blocked": len(scoped) - len(included),
        },
        "method": {
            "name": "CarnivalFlow",
            "single_joint_problem": True,
            "condition_count": len(included),
            "lambda_nominal": nominal_lambda,
            "lambda_solver": solver_lambda,
            "lambda_scaling": lambda_scaling,
            "lambda_scaling_formula": scaling_formula,
            "objective_interpretation": "mean condition fit + nominal lambda * union edge count",
            "max_flow": max_flow,
            "max_seconds": max_seconds,
        },
        "conditions": [],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        import corneto as cn
        import numpy as np
        from corneto.methods.carnival import CarnivalFlow

        graph = cn.Graph.from_tuples(
            [
                (str(edge["source"]), int(edge["sign"]), str(edge["target"]))
                for edge in graph_rows
            ]
        )
        cdict: dict[str, dict[str, dict[str, object]]] = {}
        for row in included:
            features = {
                str(node): {"value": float(value), "mapping": "vertex", "role": "input"}
                for node, value in row["inputs"].items()
            }
            features.update(
                {
                    str(node): {
                        "value": float(value),
                        "mapping": "vertex",
                        "role": "output",
                    }
                    for node, value in row["outputs"].items()
                }
            )
            cdict[str(row["run_accession"])] = features
        data = cn.Data.from_cdict(cdict)
        method = CarnivalFlow(
            lambda_reg=solver_lambda,
            max_flow=max_flow,
            enable_bfs_heuristic=True,
        )
        problem = method.build_from_data(graph, data)
        processed_names = list(method.processed_data.samples)
        if processed_names != condition_names:
            raise ValueError(
                f"CORNETO condition order changed: {processed_names!r} != {condition_names!r}"
            )
        selected_solver = solver
        if selected_solver == "auto":
            available = [str(value).upper() for value in cn.opt.available_solvers()]
            selected_solver = "gurobi" if "GUROBI" in available else "highs"
        solved = problem.solve(
            solver=selected_solver,
            max_seconds=max_seconds,
            verbosity=0,
        )
        solver_status = str(getattr(solved, "status", "unknown"))
        edge_expr = problem.expr.get("edge_has_signal")
        selected_by_condition: dict[str, list[dict[str, object]]] = {
            name: [] for name in condition_names
        }
        has_incumbent = edge_expr is not None and edge_expr.value is not None
        if has_incumbent:
            selected_by_condition = _edge_records(method, edge_expr.value, condition_names)

        vertex_expr = problem.expr.get("vertex_value")
        vertex_values = None
        if vertex_expr is not None and vertex_expr.value is not None:
            vertex_values = np.asarray(vertex_expr.value, dtype=float)
            if vertex_values.ndim == 1:
                vertex_values = vertex_values.reshape((-1, 1))
            if (
                vertex_values.shape[1] != len(condition_names)
                and vertex_values.shape[0] == len(condition_names)
            ):
                vertex_values = vertex_values.T
        condition_receipts: list[dict[str, object]] = []
        for row in scoped:
            run = str(row["run_accession"])
            if row.get("preprocessing_status") != "included":
                condition_receipts.append(
                    {
                        "study_accession": row.get("study_accession"),
                        "run_accession": run,
                        "canonical_ocm_id": row.get("canonical_ocm_id"),
                        "patient_id": row.get("patient_id"),
                        "status": "preprocessing_blocked",
                        "reason": row.get("blocked_reason"),
                        "selected_edges": [],
                    }
                )
                continue
            condition_index = condition_names.index(run)
            agreements: list[bool] = []
            if vertex_values is not None:
                for node, expected in row["outputs"].items():
                    if node not in method.processed_graph.V:
                        continue
                    vertex_index = method.processed_graph.V.index(node)
                    observed = float(vertex_values[vertex_index, condition_index])
                    agreements.append(observed * float(expected) > 0)
            condition_receipts.append(
                {
                    "study_accession": row.get("study_accession"),
                    "run_accession": run,
                    "canonical_ocm_id": row.get("canonical_ocm_id"),
                    "patient_id": row.get("patient_id"),
                    "status": solver_status,
                    "selected_edges": selected_by_condition[run],
                    "selected_edge_count": len(selected_by_condition[run]),
                    "output_sign_agreement": (
                        sum(agreements) / len(agreements) if agreements else None
                    ),
                    "output_count_evaluated": len(agreements),
                }
            )
        union_edges = {
            (str(edge["source"]), str(edge["target"]), int(edge.get("sign", 0)))
            for edges in selected_by_condition.values()
            for edge in edges
        }
        receipt.update(
            {
                "status": "completed" if has_incumbent else "no_incumbent",
                "solver": {
                    "requested": solver,
                    "selected": selected_solver.upper(),
                    "status": solver_status,
                    "has_incumbent": has_incumbent,
                    "objective_value": _json_value(getattr(solved, "value", None)),
                },
                "software": {
                    "python": platform.python_version(),
                    "corneto": getattr(cn, "__version__", "unknown"),
                },
                "processed_graph": {
                    "vertices": len(method.processed_graph.V),
                    "edges": len(method.processed_graph.E),
                },
                "selected_edge_union_count": len(union_edges),
                "conditions": condition_receipts,
                "claim_limit": (
                    "response-blind joint network stability only; no drug-response or causal "
                    "interpretation"
                ),
            }
        )
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error_class"] = type(exc).__name__
        receipt["error_message"] = str(exc)[:1000]
        raise
    finally:
        receipt["finished_at_utc"] = datetime.now(UTC).replace(microsecond=0).isoformat()
        output.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "mode": mode,
                "study": study,
                "lambda_nominal": nominal_lambda,
                "lambda_solver": solver_lambda,
                "conditions": len(included),
                "output": str(output),
            },
            sort_keys=True,
        )
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("pooled", "cohort"), required=True)
    parser.add_argument("--study")
    parser.add_argument("--lambda-reg", type=float, required=True, dest="nominal_lambda")
    parser.add_argument("--lambda-scaling", choices=("mean_fit", "raw"), default="mean_fit")
    parser.add_argument("--solver", choices=("auto", "gurobi", "highs"), default="auto")
    parser.add_argument("--max-flow", type=int, default=100)
    parser.add_argument("--max-seconds", type=int, default=259200)
    parser.add_argument("--max-conditions", type=int)
    args = parser.parse_args()
    return solve_bundle(
        bundle_path=args.bundle,
        output=args.output,
        mode=args.mode,
        study=args.study,
        nominal_lambda=args.nominal_lambda,
        lambda_scaling=args.lambda_scaling,
        solver=args.solver,
        max_flow=args.max_flow,
        max_seconds=args.max_seconds,
        max_conditions=args.max_conditions,
    )


if __name__ == "__main__":
    raise SystemExit(main())
