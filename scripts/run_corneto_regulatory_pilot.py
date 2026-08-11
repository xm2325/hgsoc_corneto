#!/usr/bin/env python3
"""Response-blind TF-activity/CARNIVAL pilot using public regulatory priors.

This is deliberately an auditable pilot, not a Taxol mechanism claim.  It
computes signed regulon scores from CollecTRI, chooses deterministic extreme
TF outputs and expression-derived upstream priors, then solves CORNETO's
CarnivalFlow on a bounded OmniPath subgraph.  Missing paths or solver failures
are recorded as blocked/error; they are never silently converted to a result.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path


def _bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _find_column(fieldnames: list[str], candidates: tuple[str, ...]) -> str:
    lower = {x.lower(): x for x in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise ValueError(f"none of {candidates!r} found in columns {fieldnames!r}")


def _load_manifest(path: Path, study: str, primary_only: bool) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        run_col = _find_column(fields, ("run_accession", "run", "sample_id"))
        study_col = next((x for x in fields if x.lower() in {"study_accession", "study"}), None)
        rows = []
        for row in reader:
            if study_col and row.get(study_col, "").strip() != study:
                continue
            if primary_only:
                eligible = row.get("primary_cohort_eligible", "")
                sample_class = row.get("sample_class", row.get("sample_type", ""))
                if eligible and not _bool(eligible):
                    continue
                if sample_class and sample_class.strip().lower() not in {"tumour", "tumor", "tumour_primary", "primary_tumour"}:
                    continue
            run = row.get(run_col, "").strip()
            if run:
                rows.append({"run_accession": run, **{k: str(v or "") for k, v in row.items()}})
    dedup = {row["run_accession"]: row for row in rows}
    return [dedup[key] for key in sorted(dedup)]


def _load_expression(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            raise ValueError("expression matrix has no sample columns")
        sample_ids = header[1:]
        values: dict[str, list[float]] = {}
        for row in reader:
            if len(row) < len(header):
                continue
            gene = row[0].strip()
            if not gene or gene in values:
                continue
            out = []
            for item in row[1:]:
                try:
                    value = float(item)
                except ValueError:
                    value = 0.0
                out.append(math.log1p(max(0.0, value)))
            values[gene] = out
    return sample_ids, values


def _load_edges(path: Path) -> list[tuple[str, str, int]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    edges = set()
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            source = (row.get("source") or "").strip()
            target = (row.get("target") or "").strip()
            try:
                sign = int(float(row.get("sign", "0")))
            except ValueError:
                continue
            if source and target and sign in {-1, 1} and source != target:
                edges.add((source, target, sign))
    return sorted(edges)


def _z_scores(values: dict[str, list[float]], sample_ids: list[str]) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {}
    for gene, row in values.items():
        mean = statistics.fmean(row) if row else 0.0
        sd = statistics.pstdev(row) if len(row) > 1 else 0.0
        result[gene] = [0.0 if sd == 0 else (x - mean) / sd for x in row]
    return result


def _regulon_scores(edges, zscores, sample_index, min_targets):
    targets: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for source, target, sign in edges:
        if target in zscores:
            targets[source].append((target, sign))
    scores = {}
    for tf, tf_edges in targets.items():
        if len(tf_edges) < min_targets:
            continue
        nums = [sign * zscores[target][sample_index] for target, sign in tf_edges]
        scores[tf] = sum(nums) / math.sqrt(len(nums))
    return scores, {tf: len(v) for tf, v in targets.items() if len(v) >= min_targets}


def _select_network(all_edges, scores, zscores, sample_index, max_outputs, max_inputs, max_depth, max_edges):
    outputs = [(node, value) for node, value in scores.items() if math.isfinite(value) and abs(value) > 1.0]
    outputs.sort(key=lambda item: (-abs(item[1]), item[0]))
    outputs = outputs[:max_outputs]
    output_nodes = {node for node, _ in outputs}
    incoming = defaultdict(list)
    outgoing = defaultdict(list)
    indegree = defaultdict(int)
    nodes = set()
    for source, target, sign in all_edges:
        incoming[target].append((source, sign))
        outgoing[source].append((target, sign))
        indegree[target] += 1
        nodes.update((source, target))
    upstream = set(output_nodes)
    frontier = set(output_nodes)
    for _ in range(max_depth):
        nxt = set()
        for node in frontier:
            nxt.update(source for source, _ in incoming[node])
        nxt -= upstream
        upstream.update(nxt)
        frontier = nxt
    source_candidates = [node for node in upstream if indegree[node] == 0 and node in zscores]
    source_candidates.sort(key=lambda node: (-abs(zscores[node][sample_index]), node))
    # When an upstream component has no explicit source, use high-scoring
    # boundary nodes as expression-derived priors and record that policy.
    if not source_candidates:
        source_candidates = [node for node in upstream if node in zscores]
        source_candidates.sort(key=lambda node: (-abs(zscores[node][sample_index]), node))
    inputs = source_candidates[:max_inputs]
    relevant = set(inputs) | output_nodes
    # Keep a bounded union of paths: reverse reachability from outputs and
    # forward reachability from chosen inputs, with deterministic edge order.
    for _ in range(max_depth):
        for source, target, _ in all_edges:
            if target in relevant and source in upstream:
                relevant.add(source)
            if source in relevant and target in upstream:
                relevant.add(target)
    selected = [(s, t, sign) for s, t, sign in all_edges if s in relevant and t in relevant]
    selected = sorted(selected, key=lambda item: (item[0], item[1], item[2]))[:max_edges]
    return inputs, outputs, selected


def _solve_carnival(edges, inputs, outputs, solver, max_flow, max_seconds):
    import corneto as cn
    from corneto.methods.carnival import CarnivalFlow

    graph = cn.Graph.from_tuples(edges)
    method = CarnivalFlow(lambda_reg=1.0, max_flow=max_flow, enable_bfs_heuristic=True)
    problem = method.build(graph, perturbations=inputs, transcription_factors=outputs)
    solved = problem.solve(solver=solver, max_seconds=max_seconds, verbosity=0)
    status = str(getattr(solved, "status", "unknown"))
    edge_values = problem.expr.get("edge_has_signal")
    selected = []
    if edge_values is not None and edge_values.value is not None:
        values = edge_values.value
        for index, value in enumerate(values):
            scalar = float(value[0] if hasattr(value, "__len__") else value)
            if scalar > 0.5 and index < len(method.processed_graph.E):
                source, target = method.processed_graph.E[index]
                source = sorted(source)
                target = sorted(target)
                if source and target:
                    selected.append({"source": source[0], "target": target[0], "sign": method.processed_graph.get_attr_edge(index).get("interaction")})
    return {"status": status, "objective": getattr(solved, "value", None),
            "selected_edges": selected, "processed_vertices": len(method.processed_graph.V),
            "processed_edges": len(method.processed_graph.E)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--collectri", type=Path, required=True)
    parser.add_argument("--pkn", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primary-only", action="store_true")
    parser.add_argument("--max-samples", type=int, default=8)
    parser.add_argument("--min-targets", type=int, default=5)
    parser.add_argument("--max-inputs", type=int, default=5)
    parser.add_argument("--max-outputs", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=3)
    parser.add_argument("--max-edges", type=int, default=5000)
    parser.add_argument("--max-flow", type=int, default=20)
    parser.add_argument("--max-seconds", type=int, default=600)
    parser.add_argument("--solver", choices=("auto", "gurobi", "highs"), default="auto")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema_version": "corneto_regulatory_pilot.v1",
        "status": "running",
        "started_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "study": args.study, "primary_only": args.primary_only,
        "method": {"name": "CarnivalFlow", "response_blind": True,
                    "activity_method": "signed_collecTRI_regulon_zscore",
                    "input_policy": "expression-derived upstream prior; not a perturbation experiment",
                    "min_targets": args.min_targets, "max_depth": args.max_depth,
                    "max_inputs": args.max_inputs, "max_outputs": args.max_outputs,
                    "max_edges": args.max_edges, "max_flow": args.max_flow},
        "inputs": {"expression": str(args.expression), "manifest": str(args.manifest),
                   "collectri": str(args.collectri), "pkn": str(args.pkn)},
        "source_sha256": {key: _sha(path) for key, path in (("expression", args.expression), ("manifest", args.manifest), ("collectri", args.collectri), ("pkn", args.pkn))},
        "samples": [],
    }
    try:
        manifest_rows = _load_manifest(args.manifest, args.study, args.primary_only)
        expression_samples, values = _load_expression(args.expression)
        zscores = _z_scores(values, expression_samples)
        edges_collectri = _load_edges(args.collectri)
        edges_pkn = _load_edges(args.pkn)
        selected_manifest = [row for row in manifest_rows if row["run_accession"] in expression_samples][:args.max_samples]
        if not selected_manifest:
            raise ValueError("no manifest samples overlap expression columns")
        receipt["input_counts"] = {"manifest_rows": len(manifest_rows), "expression_samples": len(expression_samples),
                                    "expression_genes": len(values), "collectri_edges": len(edges_collectri), "pkn_edges": len(edges_pkn),
                                    "selected_samples": len(selected_manifest)}
        solver = args.solver
        if solver == "auto":
            try:
                import corneto as cn
                solver = "gurobi" if "GUROBI" in [x.upper() for x in cn.opt.available_solvers()] else "highs"
            except Exception:
                solver = "highs"
        receipt["solver"] = {"requested": args.solver, "selected": solver.upper()}
        sample_index = {sample: index for index, sample in enumerate(expression_samples)}
        tf_all = defaultdict(list)
        for source, target, sign in edges_collectri:
            tf_all[source].append((target, sign))
        for row in selected_manifest:
            run = row["run_accession"]
            index = sample_index[run]
            scores, target_counts = _regulon_scores(edges_collectri, zscores, index, args.min_targets)
            # Recompute source ranking for this sample (the helper uses the first
            # column only for fallback source ranking, so replace it below).
            inputs, outputs, network = _select_network(edges_pkn, scores, zscores, index, args.max_outputs, args.max_inputs, args.max_depth, args.max_edges)
            input_values = {node: (1.0 if zscores.get(node, [0.0] * len(expression_samples))[index] >= 0 else -1.0) for node in inputs}
            output_values = {node: (1.0 if value >= 0 else -1.0) for node, value in outputs}
            sample_result = {"run_accession": run, "regulon_count": len(scores),
                             "top_tf_scores": [{"tf": tf, "score": score, "target_count": target_counts.get(tf, 0)} for tf, score in outputs],
                             "inputs": input_values, "outputs": output_values,
                             "network_edges": len(network), "status": "blocked", "selected_edges": []}
            if not network or not input_values or not output_values:
                sample_result["reason"] = "no_bounded_signed_path_between_expression_priors_and_TFs"
            else:
                try:
                    sample_result.update(_solve_carnival(network, input_values, output_values, solver, args.max_flow, args.max_seconds))
                except Exception as exc:
                    sample_result["status"] = "error"
                    sample_result["error_class"] = type(exc).__name__
                    sample_result["error_message"] = str(exc)[:500]
            receipt["samples"].append(sample_result)
        statuses = [row["status"] for row in receipt["samples"]]
        receipt["status"] = "completed" if any(status in {"optimal", "optimal_inaccurate"} for status in statuses) else "blocked"
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error_class"] = type(exc).__name__
        receipt["error_message"] = str(exc)[:500]
        raise
    finally:
        receipt["finished_at_utc"] = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
