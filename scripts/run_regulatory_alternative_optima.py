#!/usr/bin/env python3
"""Response-blind alternative-solution stability for regulatory CORNETO.

The driver rebuilds the same deterministic, per-sample CARNIVAL problems used
by ``run_corneto_regulatory_pilot.py`` and samples near-optimal solutions by
perturbing ``edge_has_signal``.  It never reads a drug-response phenotype.

The CORNETO sampler always returns the incumbent as solution zero; subsequent
solutions are accepted perturbations whose non-excluded objectives remain
within ``rel_opt_tol`` of the incumbent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

from run_corneto_regulatory_pilot import (
    _load_edges,
    _load_expression,
    _load_manifest,
    _regulon_scores,
    _select_network,
    _z_scores,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _edge_id(source: str, target: str, sign: Any) -> str:
    return f"{source}|{target}|{sign}"


def _active_edge_sets(values: np.ndarray, processed_graph: Any) -> list[set[str]]:
    """Convert sampler output (solutions, edges, ...) to signed edge sets."""
    array = np.asarray(values)
    if array.ndim < 2:
        raise ValueError(f"edge_has_signal samples must be at least 2D, got {array.shape}")
    edge_count = len(processed_graph.E)
    if array.shape[1] != edge_count:
        raise ValueError(
            f"edge_has_signal edge axis={array.shape[1]} but processed graph has {edge_count} edges"
        )
    result: list[set[str]] = []
    for solution in array:
        active = np.asarray(solution).reshape(edge_count, -1).max(axis=1) > 0.5
        edges: set[str] = set()
        for index in np.flatnonzero(active):
            source_set, target_set = processed_graph.E[int(index)]
            sources = sorted(source_set)
            targets = sorted(target_set)
            if not sources or not targets:
                continue
            sign = processed_graph.get_attr_edge(int(index)).get("interaction")
            edges.add(_edge_id(sources[0], targets[0], sign))
        result.append(edges)
    return result


def _jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def summarize_edge_ensemble(edge_sets: list[set[str]], core_frequency: float) -> dict[str, Any]:
    """Return deterministic ensemble metrics; solution zero is the incumbent."""
    if not edge_sets:
        raise ValueError("edge ensemble is empty")
    if not 0.0 < core_frequency <= 1.0:
        raise ValueError("core_frequency must be in (0, 1]")
    solution_count = len(edge_sets)
    union = set().union(*edge_sets)
    intersection = set.intersection(*edge_sets)
    frequencies = {
        edge: sum(edge in solution for solution in edge_sets) / solution_count
        for edge in sorted(union)
    }
    pairwise = [_jaccard(a, b) for a, b in combinations(edge_sets, 2)]
    incumbent_jaccard = [_jaccard(edge_sets[0], other) for other in edge_sets[1:]]
    entropies = []
    for frequency in frequencies.values():
        entropy = 0.0
        for probability in (frequency, 1.0 - frequency):
            if probability > 0:
                entropy -= probability * math.log2(probability)
        entropies.append(entropy)
    return {
        "solution_count": solution_count,
        "accepted_alternative_count": solution_count - 1,
        "incumbent_edge_count": len(edge_sets[0]),
        "edge_union_count": len(union),
        "edge_intersection_count": len(intersection),
        "core_edge_count": sum(value >= core_frequency for value in frequencies.values()),
        "core_frequency_threshold": core_frequency,
        "mean_pairwise_jaccard": None if not pairwise else sum(pairwise) / len(pairwise),
        "min_pairwise_jaccard": None if not pairwise else min(pairwise),
        "mean_incumbent_jaccard": (
            None if not incumbent_jaccard else sum(incumbent_jaccard) / len(incumbent_jaccard)
        ),
        "mean_edge_entropy_bits": None if not entropies else sum(entropies) / len(entropies),
        "edge_frequencies": [
            {"edge": edge, "frequency": frequency}
            for edge, frequency in sorted(frequencies.items(), key=lambda item: (-item[1], item[0]))
        ],
    }


def _build_problem(edges, inputs, outputs, lambda_reg, max_flow):
    import corneto as cn
    from corneto.methods.carnival import CarnivalFlow

    graph = cn.Graph.from_tuples([(source, sign, target) for source, target, sign in edges])
    method = CarnivalFlow(
        lambda_reg=lambda_reg,
        max_flow=max_flow,
        enable_bfs_heuristic=True,
    )
    features = {
        node: {"value": value, "mapping": "vertex", "role": "input"}
        for node, value in inputs.items()
    }
    features.update(
        {
            node: {"value": value, "mapping": "vertex", "role": "output"}
            for node, value in outputs.items()
        }
    )
    data = cn.Data.from_cdict({"condition": features})
    return method, method.build_from_data(graph, data)


def _objective_snapshot(problem) -> list[dict[str, Any]]:
    snapshot = []
    for objective in problem.objectives:
        value = objective.value
        try:
            scalar = float(np.asarray(value).reshape(-1)[0])
        except (TypeError, ValueError, IndexError):
            scalar = None
        snapshot.append({"name": str(objective.name), "value": scalar})
    return snapshot


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--collectri", type=Path, required=True)
    parser.add_argument("--pkn", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primary-only", action="store_true")
    parser.add_argument("--max-samples", type=int, default=1)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--sample-index",
        type=int,
        help="Zero-based index in the deterministic study/run-sorted eligible sample list.",
    )
    selection.add_argument(
        "--sample-run",
        help="Run accession to select from the eligible expression/manifest overlap.",
    )
    parser.add_argument("--min-targets", type=int, default=5)
    parser.add_argument("--max-inputs", type=int, default=5)
    parser.add_argument("--max-outputs", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=3)
    parser.add_argument("--max-edges", type=int, default=5000)
    parser.add_argument("--max-flow", type=int, default=20)
    parser.add_argument("--lambda-reg", type=float, default=0.0)
    parser.add_argument("--percentage", type=float, default=0.03)
    parser.add_argument("--scale", type=float, default=0.03)
    parser.add_argument("--rel-opt-tol", type=float, default=0.05)
    parser.add_argument("--max-alternatives", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--core-frequency", type=float, default=0.8)
    parser.add_argument("--max-seconds", type=int, default=300)
    parser.add_argument("--solver", choices=("gurobi", "highs"), default="gurobi")
    parser.add_argument(
        "--exclude-objectives-pattern",
        default="regularization",
        help="Regex passed to CORNETO sampler; use NONE to check every objective.",
    )
    args = parser.parse_args()
    if args.max_samples < 1 or args.max_alternatives < 1:
        parser.error("max-samples and max-alternatives must be positive")
    if not 0 < args.percentage <= 1:
        parser.error("percentage must be in (0, 1]")
    if args.scale <= 0 or args.rel_opt_tol < 0:
        parser.error("scale must be positive and rel-opt-tol non-negative")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    receipt: dict[str, Any] = {
        "schema_version": "regulatory_alternative_optima.v1",
        "status": "running",
        "started_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "study": args.study,
        "response_blind": True,
        "phenotype_inputs": [],
        "method": {
            "network_method": "CarnivalFlow",
            "sampler": "corneto.methods.sampler.sample_alternative_solutions",
            "perturbed_variable": "edge_has_signal",
            "lambda_reg": args.lambda_reg,
            "percentage": args.percentage,
            "scale": args.scale,
            "rel_opt_tol": args.rel_opt_tol,
            "max_alternative_trials": args.max_alternatives,
            "seed": args.seed,
            "exclude_objectives_pattern": (
                None
                if args.exclude_objectives_pattern.upper() == "NONE"
                else args.exclude_objectives_pattern
            ),
        },
        "inputs": {
            "expression": str(args.expression),
            "manifest": str(args.manifest),
            "collectri": str(args.collectri),
            "pkn": str(args.pkn),
        },
        "source_sha256": {
            label: _sha256(path)
            for label, path in (
                ("expression", args.expression),
                ("manifest", args.manifest),
                ("collectri", args.collectri),
                ("pkn", args.pkn),
            )
        },
        "samples": [],
    }
    try:
        from corneto.methods.sampler import sample_alternative_solutions

        manifest_rows = _load_manifest(args.manifest, args.study, args.primary_only)
        expression_samples, expression = _load_expression(args.expression)
        zscores = _z_scores(expression, expression_samples)
        collectri = _load_edges(args.collectri)
        pkn = _load_edges(args.pkn)
        eligible = [
            row for row in manifest_rows if row["run_accession"] in set(expression_samples)
        ]
        if args.sample_index is not None:
            if args.sample_index < 0 or args.sample_index >= len(eligible):
                raise ValueError(
                    f"sample-index {args.sample_index} outside eligible range 0..{len(eligible) - 1}"
                )
            selected = [eligible[args.sample_index]]
        elif args.sample_run is not None:
            selected = [row for row in eligible if row["run_accession"] == args.sample_run]
            if not selected:
                raise ValueError(f"sample-run {args.sample_run!r} is not eligible for {args.study}")
        else:
            selected = eligible[: args.max_samples]
        if not selected:
            raise ValueError("no manifest samples overlap expression columns")
        sample_index = {sample: index for index, sample in enumerate(expression_samples)}
        receipt["input_counts"] = {
            "manifest_rows": len(manifest_rows),
            "expression_samples": len(expression_samples),
            "expression_genes": len(expression),
            "collectri_edges": len(collectri),
            "pkn_edges": len(pkn),
            "eligible_samples": len(eligible),
            "selected_samples": len(selected),
        }
        for sample_number, row in enumerate(selected):
            run = row["run_accession"]
            index = sample_index[run]
            scores, target_counts = _regulon_scores(collectri, zscores, index, args.min_targets)
            inputs, outputs, network = _select_network(
                pkn,
                scores,
                zscores,
                index,
                args.max_outputs,
                args.max_inputs,
                args.max_depth,
                args.max_edges,
            )
            input_values = {
                node: 1.0 if zscores[node][index] >= 0 else -1.0 for node in inputs
            }
            output_values = {node: 1.0 if value >= 0 else -1.0 for node, value in outputs}
            present = {node for edge in network for node in edge[:2]}
            input_values = {node: value for node, value in input_values.items() if node in present}
            output_values = {node: value for node, value in output_values.items() if node in present}
            result: dict[str, Any] = {
                "run_accession": run,
                "network_edges": len(network),
                "inputs": input_values,
                "outputs": output_values,
                "top_tf_scores": [
                    {"tf": tf, "score": score, "target_count": target_counts.get(tf, 0)}
                    for tf, score in outputs
                ],
                "status": "blocked",
            }
            if not network or not input_values or not output_values:
                result["reason"] = "no_bounded_signed_path_between_expression_priors_and_TFs"
                receipt["samples"].append(result)
                continue
            try:
                method, problem = _build_problem(
                    network, input_values, output_values, args.lambda_reg, args.max_flow
                )
                baseline = problem.solve(
                    solver=args.solver,
                    max_seconds=args.max_seconds,
                    verbosity=0,
                )
                baseline_status = str(getattr(baseline, "status", "unknown"))
                result["baseline_status"] = baseline_status
                result["baseline_objectives"] = _objective_snapshot(problem)
                if baseline_status not in {"optimal", "optimal_inaccurate"}:
                    result["status"] = "blocked_baseline_not_optimal"
                else:
                    sampled = sample_alternative_solutions(
                        problem,
                        "edge_has_signal",
                        percentage=args.percentage,
                        scale=args.scale,
                        rel_opt_tol=args.rel_opt_tol,
                        max_samples=args.max_alternatives,
                        solver_kwargs={
                            "solver": args.solver,
                            "max_seconds": args.max_seconds,
                        },
                        rng=args.seed + sample_number,
                        collect_vars=["edge_has_signal"],
                        exclude_objectives_pattern=(
                            None
                            if args.exclude_objectives_pattern.upper() == "NONE"
                            else args.exclude_objectives_pattern
                        ),
                        verbose=0,
                    )
                    edge_sets = _active_edge_sets(sampled["edge_has_signal"], method.processed_graph)
                    result.update(summarize_edge_ensemble(edge_sets, args.core_frequency))
                    result["status"] = (
                        "completed" if edge_sets[0] else "blocked_no_selected_edges"
                    )
                    result["processed_vertices"] = len(method.processed_graph.V)
                    result["processed_edges"] = len(method.processed_graph.E)
            except Exception as error:  # fail closed but retain per-sample evidence
                result["status"] = "error"
                result["error_class"] = type(error).__name__
                result["error_message"] = str(error)[:1000]
            receipt["samples"].append(result)
        statuses = [sample["status"] for sample in receipt["samples"]]
        receipt["status_counts"] = {
            status: statuses.count(status) for status in sorted(set(statuses))
        }
        receipt["status"] = "completed" if "completed" in statuses else "blocked"
    except Exception as error:
        receipt["status"] = "failed"
        receipt["error_class"] = type(error).__name__
        receipt["error_message"] = str(error)[:1000]
        raise
    finally:
        receipt["finished_at_utc"] = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        args.output.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(receipt, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
