#!/usr/bin/env python3
"""Response-blind deterministic regulatory pilot.

This pilot computes signed CollecTRI regulon scores and a bounded, normalized
OmniPath signed subgraph. It deliberately does not claim a CORNETO/CARNIVAL
solution: until a pinned adapter/API is verified, the receipt records an
explicit blocked_unverified_api status instead of fabricating solver output.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def _bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalise_symbol(value: str | None) -> str | None:
    text = (value or "").strip().upper()
    if not text or text in {"NA", "N/A", "NULL", "NONE"}:
        return None
    if any(char in text for char in (":", ";", "|")):
        return None
    return text


def _load_manifest(path: Path, study: str, primary_only: bool) -> list[dict[str, str]]:
    required = {"run_accession", "study_accession"}
    if primary_only:
        required |= {
            "sample_class",
            "histotype_group",
            "is_representative_rna_library",
            "primary_cohort_eligible",
        }
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(f"manifest is missing required columns: {missing}")
        rows: list[dict[str, str]] = []
        for source_row in reader:
            row = {key: str(value or "") for key, value in source_row.items()}
            if row.get("study_accession", "").strip() != study:
                continue
            if primary_only:
                sample_class = row.get("sample_class", "").strip().lower()
                if sample_class not in {"tumour", "tumor"}:
                    continue
                if row.get("histotype_group", "").strip().upper() != "HGSOC":
                    continue
                if not _bool(row.get("is_representative_rna_library")):
                    continue
                if not _bool(row.get("primary_cohort_eligible")):
                    continue
            run = row.get("run_accession", "").strip()
            if run:
                rows.append({"run_accession": run, **row})
    dedup = {row["run_accession"]: row for row in rows}
    return [dedup[key] for key in sorted(dedup)]


def _load_expression(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    """Read a Salmon-style gene matrix with robust gene-name detection.

    The usual header is gene_id, gene_name, sample.... Some exports use
    gene_symbol or symbol instead; those are treated identically. For a matrix
    without a name column, the historical fallback is the first column as gene
    ID and values beginning at column 1.
    """
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"empty expression matrix: {path}") from exc
        lowered = [str(value).strip().lower() for value in header]
        name_index = next(
            (index for index, value in enumerate(lowered)
             if value in {"gene_name", "gene_symbol", "symbol"}),
            None,
        )
        if name_index is None:
            sample_start = 1
            if len(header) < 2:
                raise ValueError("expression matrix needs a gene column and samples")
        else:
            # For the canonical gene_id/gene_name layout this is exactly 2.
            sample_start = name_index + 1
            if len(header) <= sample_start:
                raise ValueError("expression matrix has no sample columns")
        sample_ids = [str(value).strip() for value in header[sample_start:]]
        if any(not value for value in sample_ids):
            raise ValueError("expression matrix contains an empty sample column")
        values: dict[str, list[float]] = {}
        counts: dict[str, int] = defaultdict(int)
        for row in reader:
            if len(row) != len(header):
                continue
            raw_gene = row[name_index] if name_index is not None else row[0]
            gene = _normalise_symbol(raw_gene)
            if not gene:
                continue
            numbers: list[float] = []
            for item in row[sample_start:]:
                try:
                    value = float(item)
                except (TypeError, ValueError):
                    value = 0.0
                numbers.append(math.log1p(max(0.0, value)))
            if gene not in values:
                values[gene] = [0.0] * len(sample_ids)
            values[gene] = [
                left + right for left, right in zip(values[gene], numbers)
            ]
            counts[gene] += 1
        for gene, numbers in values.items():
            values[gene] = [number / counts[gene] for number in numbers]
    return sample_ids, values


def _edge_sign(row: dict[str, str], fields: list[str]) -> int | None:
    for name in ("consensus_direction", "direction", "sign"):
        field = next((item for item in fields if item.lower() == name), None)
        if field:
            raw = str(row.get(field, "")).strip().lower()
            try:
                value = float(raw)
            except ValueError:
                value = 0.0
            if value > 0 or raw in {"positive", "activation", "activating", "stimulation"}:
                return 1
            if value < 0 or raw in {"negative", "inhibition", "inhibiting", "inhibitory"}:
                return -1
    stimulation = next(
        (item for item in fields
         if item.lower() in {"is_stimulation", "consensus_stimulation", "stimulation"}),
        None,
    )
    inhibition = next(
        (item for item in fields
         if item.lower() in {"is_inhibition", "consensus_inhibition", "inhibition"}),
        None,
    )
    stim = _bool(row.get(stimulation, "")) if stimulation else False
    inhib = _bool(row.get(inhibition, "")) if inhibition else False
    if stim and not inhib:
        return 1
    if inhib and not stim:
        return -1
    return None


def _load_edges(path: Path) -> list[tuple[str, str, int]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    edges: set[tuple[str, str, int]] = set()
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        source_field = next(
            (item for item in fields if item.lower() in {
                "source_genesymbol", "source_gene", "source",
            }),
            None,
        )
        target_field = next(
            (item for item in fields if item.lower() in {
                "target_genesymbol", "target_gene", "target",
            }),
            None,
        )
        if source_field is None or target_field is None:
            raise ValueError(f"{path} lacks source/target columns")
        for row in reader:
            source = _normalise_symbol(row.get(source_field))
            target = _normalise_symbol(row.get(target_field))
            sign = _edge_sign(row, fields)
            if source and target and sign in {-1, 1} and source != target:
                edges.add((source, target, sign))
    if not edges:
        raise ValueError(f"no signed edges could be normalized from {path}")
    return sorted(edges)


def _z_scores(values: dict[str, list[float]]) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {}
    for gene, row in values.items():
        mean = statistics.fmean(row) if row else 0.0
        sd = statistics.pstdev(row) if len(row) > 1 else 0.0
        result[gene] = [0.0 if sd == 0 else (x - mean) / sd for x in row]
    return result


def _regulon_scores(
    edges: list[tuple[str, str, int]],
    zscores: dict[str, list[float]],
    sample_index: int,
    min_targets: int,
) -> tuple[dict[str, float], dict[str, int]]:
    targets: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for source, target, sign in edges:
        if target in zscores:
            targets[source].append((target, sign))
    scores: dict[str, float] = {}
    counts: dict[str, int] = {}
    for tf, tf_edges in targets.items():
        if len(tf_edges) < min_targets:
            continue
        numbers = [sign * zscores[target][sample_index] for target, sign in tf_edges]
        scores[tf] = sum(numbers) / math.sqrt(len(numbers))
        counts[tf] = len(tf_edges)
    return scores, counts


def _select_network(
    edges: list[tuple[str, str, int]],
    scores: dict[str, float],
    zscores: dict[str, list[float]],
    sample_index: int,
    max_outputs: int,
    max_inputs: int,
    max_depth: int,
    max_edges: int,
) -> tuple[list[str], list[tuple[str, float]], list[tuple[str, str, int]]]:
    outputs = [
        (node, value) for node, value in scores.items()
        if math.isfinite(value) and abs(value) > 1.0
    ]
    outputs.sort(key=lambda item: (-abs(item[1]), item[0]))
    outputs = outputs[:max_outputs]
    output_nodes = {node for node, _ in outputs}
    incoming: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for source, target, sign in edges:
        incoming[target].append((source, sign))
    upstream = set(output_nodes)
    frontier = set(output_nodes)
    for _ in range(max_depth):
        nxt: set[str] = set()
        for node in frontier:
            nxt.update(source for source, _ in incoming[node])
        nxt -= upstream
        upstream.update(nxt)
        frontier = nxt
    source_candidates = [node for node in upstream if node in zscores]
    source_candidates.sort(
        key=lambda node: (-abs(zscores[node][sample_index]), node)
    )
    inputs = source_candidates[:max_inputs]
    relevant = set(inputs) | output_nodes
    for _ in range(max_depth):
        for source, target, _ in edges:
            if target in relevant and source in upstream:
                relevant.add(source)
            if source in relevant and target in upstream:
                relevant.add(target)
    selected = [
        (source, target, sign) for source, target, sign in edges
        if source in relevant and target in relevant
    ]
    return inputs, outputs, sorted(selected)[:max_edges]


def _normalize_network(
    network: list[tuple[str, str, int]],
    zscores: dict[str, list[float]],
    sample_index: int,
    inputs: list[str],
) -> tuple[list[dict[str, object]], float]:
    scale = max(
        (abs(zscores[source][sample_index]) for source in inputs if source in zscores),
        default=0.0,
    ) or 1.0
    normalized: list[dict[str, object]] = []
    for source, target, sign in network:
        source_score = (
            zscores[source][sample_index] if source in zscores else 0.0
        )
        normalized.append({
            "source": source,
            "target": target,
            "sign": sign,
            "source_zscore": source_score,
            "normalized_edge_weight": sign * source_score / scale,
        })
    return normalized, scale


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--collectri", type=Path, required=True)
    parser.add_argument("--pkn", "--omnipath", dest="pkn", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--primary-only", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument("--max-samples", type=int, default=8)
    parser.add_argument("--min-targets", type=int, default=5)
    parser.add_argument("--max-inputs", type=int, default=5)
    parser.add_argument("--max-outputs", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=3)
    parser.add_argument("--max-edges", type=int, default=5000)
    # Compatibility options retained by the Roihu sbatch wrapper; the
    # deterministic pilot records the policy but does not invoke a solver.
    parser.add_argument("--top-tfs", type=int, default=50)
    parser.add_argument("--expression-transform", default="log1p_tpm")
    parser.add_argument("--solver", choices=("auto", "gurobi", "highs"), default="auto")
    args = parser.parse_args()
    if args.output is None:
        if args.output_dir is None:
            parser.error("one of --output or --output-dir is required")
        args.output = args.output_dir / "regulatory_pilot_summary.json"
    if not 8 <= args.max_samples <= 12:
        parser.error("--max-samples must be between 8 and 12")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    receipt: dict[str, object] = {
        "schema_version": "corneto_regulatory_pilot.v1",
        "status": "running",
        "started_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "study": args.study,
        "primary_only": args.primary_only,
        "method": {
            "name": "deterministic_signed_regulon_pkn_normalization",
            "response_blind": True,
            "activity_method": "mean signed target z-score",
            "input_policy": "expression-derived upstream prior; not a perturbation experiment",
            "min_targets": args.min_targets,
            "max_depth": args.max_depth,
            "max_inputs": args.max_inputs,
            "max_outputs": args.max_outputs,
            "max_edges": args.max_edges,
        },
        "inputs": {
            "expression": str(args.expression),
            "manifest": str(args.manifest),
            "collectri": str(args.collectri),
            "pkn": str(args.pkn),
        },
        "source_sha256": {
            key: _sha(path)
            for key, path in (
                ("expression", args.expression),
                ("manifest", args.manifest),
                ("collectri", args.collectri),
                ("pkn", args.pkn),
            )
        },
        "samples": [],
    }
    try:
        manifest_rows = _load_manifest(args.manifest, args.study, args.primary_only)
        expression_samples, values = _load_expression(args.expression)
        zscores = _z_scores(values)
        edges_collectri = _load_edges(args.collectri)
        edges_pkn = _load_edges(args.pkn)
        expression_set = set(expression_samples)
        selected_manifest = [
            row for row in manifest_rows
            if row["run_accession"] in expression_set
        ][:args.max_samples]
        if not selected_manifest:
            raise ValueError("no manifest samples overlap expression columns")
        receipt["input_counts"] = {
            "manifest_rows": len(manifest_rows),
            "expression_samples": len(expression_samples),
            "expression_genes": len(values),
            "collectri_edges": len(edges_collectri),
            "pkn_edges": len(edges_pkn),
            "selected_samples": len(selected_manifest),
        }
        receipt["solver"] = {
            "requested": args.solver,
            "status": "not_attempted",
            "reason": "No pinned regulatory CORNETO/CARNIVAL adapter; deterministic pilot only.",
        }
        receipt["corneto"] = {
            "status": "blocked_unverified_api",
            "reason": "A solver result is intentionally withheld until a pinned adapter is verified.",
        }
        sample_index = {sample: index for index, sample in enumerate(expression_samples)}
        for row in selected_manifest:
            run = row["run_accession"]
            index = sample_index[run]
            scores, target_counts = _regulon_scores(
                edges_collectri, zscores, index, args.min_targets
            )
            inputs, outputs, network = _select_network(
                edges_pkn, scores, zscores, index, args.max_outputs,
                args.max_inputs, args.max_depth, args.max_edges,
            )
            normalized_edges, normalization_scale = _normalize_network(
                network, zscores, index, inputs
            )
            sample_result: dict[str, object] = {
                "run_accession": run,
                "regulon_count": len(scores),
                "top_tf_scores": [
                    {
                        "tf": tf,
                        "score": score,
                        "target_count": target_counts.get(tf, 0),
                    }
                    for tf, score in outputs
                ],
                "inputs": {
                    node: (1.0 if zscores[node][index] >= 0 else -1.0)
                    for node in inputs
                },
                "outputs": {
                    node: (1.0 if score >= 0 else -1.0)
                    for node, score in outputs
                },
                "network_edges": len(network),
                "normalized_edges": normalized_edges,
                "normalization_scale": normalization_scale,
                "status": "blocked_unverified_api",
                "selected_edges": [],
                "reason": "No pinned regulatory CORNETO/CARNIVAL adapter; no solver result was produced.",
            }
            if not network or not inputs or not outputs:
                sample_result["reason"] = (
                    "no_bounded_signed_path_between_expression_priors_and_TFs; "
                    "CORNETO adapter remains blocked"
                )
            cast_samples = receipt["samples"]
            assert isinstance(cast_samples, list)
            cast_samples.append(sample_result)
        receipt["status"] = "completed_response_blind_blocked"
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error_class"] = type(exc).__name__
        receipt["error_message"] = str(exc)[:500]
        raise
    finally:
        receipt["finished_at_utc"] = datetime.now(timezone.utc).replace(
            microsecond=0
        ).isoformat()
        args.output.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(receipt, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
