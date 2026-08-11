#!/usr/bin/env python3
"""Joint multi-condition CARNIVAL over one or more HGSOC RNA cohorts.

Unlike ``run_corneto_regulatory_pilot.py``, which deliberately solves one
condition at a time, this runner builds one CORNETO ``Data`` object containing
all selected conditions and solves a single ``CarnivalFlow`` problem.  The
structured-sparsity lambda therefore acts on the union of condition-specific
edge selections, as intended by CORNETO's multi-sample formulation.

For a pooled run, expression evidence is standardized within study before the
conditions are combined.  This keeps the per-sample evidence identical between
pooled and cohort-stratified fits, so their network differences can be
attributed to joint optimization rather than to a change in normalization.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from run_corneto_regulatory_pilot import (
    _load_edges,
    _load_expression,
    _load_manifest,
    _regulon_scores,
    _select_network,
    _sha,
    _z_scores,
)


def _key_value(value: str, value_type=str) -> tuple[str, object]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"expected STUDY=VALUE, got {value!r}")
    key, raw = value.split("=", 1)
    key = key.strip()
    if not key or not raw.strip():
        raise argparse.ArgumentTypeError(f"expected STUDY=VALUE, got {value!r}")
    try:
        return key, value_type(raw.strip())
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def _edge_records(method, edge_matrix, condition_index: int) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    if edge_matrix is None:
        return records
    shape = getattr(edge_matrix, "shape", ())
    for edge_index in range(min(len(method.processed_graph.E), shape[0] if shape else 0)):
        if len(shape) == 1:
            scalar = float(edge_matrix[edge_index])
        else:
            scalar = float(edge_matrix[edge_index, condition_index])
        if scalar <= 0.5:
            continue
        source, target = method.processed_graph.E[edge_index]
        source_nodes = sorted(source)
        target_nodes = sorted(target)
        # Exclude artificial source/sink edges introduced by the flow graph.
        if not source_nodes or not target_nodes:
            continue
        records.append(
            {
                "source": source_nodes[0],
                "target": target_nodes[0],
                "sign": method.processed_graph.get_attr_edge(edge_index).get("interaction"),
            }
        )
    return records


def _graph_digest(edges: set[tuple[str, str, int]]) -> str:
    payload = "".join(f"{source}\t{target}\t{sign}\n" for source, target, sign in sorted(edges))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _write_frozen_graph(path: Path, edges: set[tuple[str, str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "source\ttarget\tsign\n"
        + "".join(f"{source}\t{target}\t{sign}\n" for source, target, sign in sorted(edges)),
        encoding="utf-8",
    )


def _read_frozen_graph(path: Path) -> set[tuple[str, str, int]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "source\ttarget\tsign":
        raise ValueError(f"invalid frozen graph header: {path}")
    edges: set[tuple[str, str, int]] = set()
    for line in lines[1:]:
        source, target, sign = line.split("\t")
        edges.add((source, target, int(sign)))
    if not edges:
        raise ValueError(f"frozen graph is empty: {path}")
    return edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--expression",
        action="append",
        required=True,
        metavar="STUDY=PATH",
        help="repeat for every cohort included in the joint fit",
    )
    parser.add_argument(
        "--expected-count",
        action="append",
        default=[],
        metavar="STUDY=N",
        help="fail closed unless the selected primary count matches N",
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--collectri", type=Path, required=True)
    parser.add_argument("--pkn", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--analysis-name", required=True)
    parser.add_argument("--fit-study", action="append", default=[])
    parser.add_argument("--frozen-graph", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--primary-only", action="store_true")
    parser.add_argument("--max-samples-per-study", type=int, default=0)
    parser.add_argument("--min-targets", type=int, default=5)
    parser.add_argument("--max-inputs", type=int, default=3)
    parser.add_argument("--max-outputs", type=int, default=6)
    parser.add_argument("--max-depth", type=int, default=3)
    parser.add_argument("--max-edges", type=int, default=3000)
    parser.add_argument("--max-flow", type=int, default=100)
    parser.add_argument(
        "--lambda-reg-reported",
        type=float,
        default=0.0,
        help="reported lambda; CORNETO receives n_fit_conditions times this value",
    )
    parser.add_argument("--max-seconds", type=int, default=259200)
    parser.add_argument("--solver", choices=("gurobi", "highs"), default="gurobi")
    args = parser.parse_args()

    expression_specs = dict(_key_value(item, Path) for item in args.expression)
    expected_counts = dict(_key_value(item, int) for item in args.expected_count)
    if len(expression_specs) != len(args.expression):
        raise SystemExit("duplicate --expression study")
    unknown_expected = set(expected_counts) - set(expression_specs)
    if unknown_expected:
        raise SystemExit(f"expected counts without expression: {sorted(unknown_expected)}")
    unknown_fit_studies = set(args.fit_study) - set(expression_specs)
    if unknown_fit_studies:
        raise SystemExit(f"fit studies without expression: {sorted(unknown_fit_studies)}")
    if args.prepare_only and args.fit_study:
        raise SystemExit("--prepare-only cannot be combined with --fit-study")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    receipt: dict[str, object] = {
        "schema_version": "corneto_regulatory_joint.v1",
        "status": "running",
        "started_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "analysis_name": args.analysis_name,
        "studies": list(expression_specs),
        "primary_only": args.primary_only,
        "method": {
            "name": "CarnivalFlow",
            "joint_multi_condition": True,
            "joint_multi_sample": True,
            "response_blind": True,
            "activity_method": "signed_collecTRI_regulon_zscore",
            "zscore_scope": "within_study_before_joint_fit",
            "input_policy": "expression-derived upstream prior; not a perturbation experiment",
            "common_graph_policy": "frozen_union_of_all60_deterministic_sample_bounded_paths",
            "graph_evidence_studies": list(expression_specs),
            "fit_studies": args.fit_study or list(expression_specs),
            "lambda_reg_reported": args.lambda_reg_reported,
            "lambda_scaling": "lambda_implemented = n_fit_conditions * lambda_reported",
            "min_targets": args.min_targets,
            "max_depth": args.max_depth,
            "max_inputs": args.max_inputs,
            "max_outputs": args.max_outputs,
            "max_edges_per_sample": args.max_edges,
            "max_flow": args.max_flow,
            "exclusive_signal_paths": True,
        },
        "inputs": {
            "expression": {key: str(path) for key, path in expression_specs.items()},
            "manifest": str(args.manifest),
            "collectri": str(args.collectri),
            "pkn": str(args.pkn),
            "frozen_graph": str(args.frozen_graph),
        },
        "source_sha256": {
            "expression": {key: _sha(path) for key, path in expression_specs.items()},
            "manifest": _sha(args.manifest),
            "collectri": _sha(args.collectri),
            "pkn": _sha(args.pkn),
        },
        "samples": [],
    }

    try:
        edges_collectri = _load_edges(args.collectri)
        edges_pkn = _load_edges(args.pkn)
        condition_features: dict[str, dict[str, dict[str, object]]] = {}
        sample_metadata: dict[str, dict[str, object]] = {}
        common_edges: set[tuple[str, str, int]] = set()
        study_counts: dict[str, int] = {}
        expression_counts: dict[str, dict[str, int]] = {}

        for study, expression_path in expression_specs.items():
            expression_samples, values = _load_expression(expression_path)
            zscores = _z_scores(values, expression_samples)
            manifest_rows = _load_manifest(args.manifest, study, args.primary_only)
            selected = [row for row in manifest_rows if row["run_accession"] in expression_samples]
            if args.max_samples_per_study > 0:
                selected = selected[: args.max_samples_per_study]
            if not selected:
                raise ValueError(f"{study}: no manifest samples overlap expression columns")
            if study in expected_counts and len(selected) != expected_counts[study]:
                raise ValueError(
                    f"{study}: selected {len(selected)} primary samples, expected {expected_counts[study]}"
                )
            study_counts[study] = len(selected)
            expression_counts[study] = {
                "manifest_rows": len(manifest_rows),
                "expression_samples": len(expression_samples),
                "expression_genes": len(values),
                "selected_samples": len(selected),
            }
            sample_index = {sample: index for index, sample in enumerate(expression_samples)}
            for row in selected:
                run = row["run_accession"]
                if run in condition_features:
                    raise ValueError(f"duplicate run accession across studies: {run}")
                index = sample_index[run]
                scores, target_counts = _regulon_scores(edges_collectri, zscores, index, args.min_targets)
                inputs, outputs, network = _select_network(
                    edges_pkn,
                    scores,
                    zscores,
                    index,
                    args.max_outputs,
                    args.max_inputs,
                    args.max_depth,
                    args.max_edges,
                )
                input_values = {
                    node: (1.0 if zscores[node][index] >= 0 else -1.0)
                    for node in inputs
                    if node in zscores
                }
                output_values = {node: (1.0 if value >= 0 else -1.0) for node, value in outputs}
                present_nodes = {node for edge in network for node in edge[:2]}
                input_values = {node: value for node, value in input_values.items() if node in present_nodes}
                output_values = {node: value for node, value in output_values.items() if node in present_nodes}
                if not network or not input_values or not output_values:
                    raise ValueError(f"{study}/{run}: no bounded signed path for joint input")
                common_edges.update(network)
                features: dict[str, dict[str, object]] = {
                    node: {"value": value, "mapping": "vertex", "role": "input"}
                    for node, value in input_values.items()
                }
                features.update(
                    {
                        node: {"value": value, "mapping": "vertex", "role": "output"}
                        for node, value in output_values.items()
                    }
                )
                condition_features[run] = features
                sample_metadata[run] = {
                    "run_accession": run,
                    "study": study,
                    "regulon_count": len(scores),
                    "candidate_network_edges": len(network),
                    "inputs": input_values,
                    "outputs": output_values,
                    "top_tf_scores": [
                        {"tf": tf, "score": score, "target_count": target_counts.get(tf, 0)}
                        for tf, score in outputs
                    ],
                }

        if len(condition_features) != sum(study_counts.values()):
            raise ValueError("condition count does not equal selected study counts")
        if not common_edges:
            raise ValueError("common graph is empty")

        derived_graph_digest = _graph_digest(common_edges)
        if args.prepare_only:
            if args.frozen_graph.exists():
                raise FileExistsError(f"refusing to replace frozen graph: {args.frozen_graph}")
            _write_frozen_graph(args.frozen_graph, common_edges)
            receipt["input_counts"] = {
                "graph_evidence_study_samples": study_counts,
                "graph_evidence_samples_total": len(condition_features),
                "per_study": expression_counts,
                "collectri_edges": len(edges_collectri),
                "pkn_edges": len(edges_pkn),
                "common_candidate_edges": len(common_edges),
            }
            receipt["frozen_graph"] = {
                "path": str(args.frozen_graph),
                "sha256": derived_graph_digest,
                "edges": len(common_edges),
            }
            receipt["status"] = "completed"
            return 0

        frozen_edges = _read_frozen_graph(args.frozen_graph)
        frozen_graph_digest = _graph_digest(frozen_edges)
        if frozen_graph_digest != derived_graph_digest:
            raise ValueError(
                "frozen graph differs from graph re-derived from the declared all-primary inputs"
            )
        common_edges = frozen_edges

        fit_studies = set(args.fit_study) if args.fit_study else set(expression_specs)
        condition_features = {
            run: features
            for run, features in condition_features.items()
            if str(sample_metadata[run]["study"]) in fit_studies
        }
        if not condition_features:
            raise ValueError("no conditions remain after --fit-study filtering")
        fit_study_counts = {
            study: sum(str(sample_metadata[run]["study"]) == study for run in condition_features)
            for study in expression_specs
            if study in fit_studies
        }
        lambda_implemented = args.lambda_reg_reported * len(condition_features)
        receipt["method"]["lambda_reg_implemented"] = lambda_implemented  # type: ignore[index]
        receipt["frozen_graph"] = {
            "path": str(args.frozen_graph),
            "sha256": frozen_graph_digest,
            "edges": len(common_edges),
            "verified_against_declared_all_primary_inputs": True,
        }

        import corneto as cn
        from corneto.methods.carnival import CarnivalFlow

        graph = cn.Graph.from_tuples(
            [(source, sign, target) for source, target, sign in sorted(common_edges)]
        )
        data = cn.Data.from_cdict(condition_features)
        method = CarnivalFlow(
            lambda_reg=lambda_implemented,
            max_flow=args.max_flow,
            enable_bfs_heuristic=True,
        )
        problem = method.build_from_data(graph, data)
        condition_order = list(method.processed_data.samples)
        if set(condition_order) != set(condition_features):
            raise ValueError("CORNETO processed condition set differs from input")
        solved = problem.solve(
            solver=args.solver,
            max_seconds=args.max_seconds,
            verbosity=0,
        )
        solve_status = str(getattr(solved, "status", "unknown"))
        edge_expr = problem.expr.get("edge_has_signal")
        edge_values = None if edge_expr is None else edge_expr.value
        for condition_index, run in enumerate(condition_order):
            selected_edges = _edge_records(method, edge_values, condition_index)
            row = dict(sample_metadata[run])
            row["selected_edges"] = selected_edges
            row["status"] = (
                solve_status
                if selected_edges or solve_status not in {"optimal", "optimal_inaccurate"}
                else "blocked_no_selected_edges"
            )
            receipt["samples"].append(row)  # type: ignore[union-attr]
        receipt["input_counts"] = {
            "graph_evidence_study_samples": study_counts,
            "fit_study_samples": fit_study_counts,
            "graph_evidence_samples_total": sum(study_counts.values()),
            "selected_samples_total": len(condition_features),
            "per_study": expression_counts,
            "collectri_edges": len(edges_collectri),
            "pkn_edges": len(edges_pkn),
            "common_candidate_edges": len(common_edges),
            "processed_vertices": len(method.processed_graph.V),
            "processed_edges": len(method.processed_graph.E),
        }
        receipt["solver"] = {
            "selected": args.solver.upper(),
            "status": solve_status,
            "objective": getattr(solved, "value", None),
        }
        receipt["status"] = (
            "completed" if solve_status in {"optimal", "optimal_inaccurate"} else "blocked"
        )
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error_class"] = type(exc).__name__
        receipt["error_message"] = str(exc)[:1000]
        raise
    finally:
        receipt["finished_at_utc"] = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        args.output.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "analysis_name": args.analysis_name,
                "samples": len(receipt["samples"]),
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
