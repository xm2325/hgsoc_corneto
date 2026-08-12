#!/usr/bin/env python3
"""Freeze one independently normalized external multi-condition CARNIVAL bundle."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

try:
    from run_corneto_regulatory_pilot import (
        _load_edges,
        _load_expression,
        _load_manifest,
        _regulon_scores,
        _select_network,
        _sha,
        _z_scores,
    )
except ModuleNotFoundError:
    from scripts.run_corneto_regulatory_pilot import (
        _load_edges,
        _load_expression,
        _load_manifest,
        _regulon_scores,
        _select_network,
        _sha,
        _z_scores,
    )


def build_bundle(args: argparse.Namespace) -> dict[str, object]:
    for path in (args.expression, args.manifest, args.collectri, args.pkn):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output.exists():
        raise ValueError(f"refusing to overwrite {args.output}")
    manifest_rows = _load_manifest(args.manifest, args.study, primary_only=False)
    include_roles = set(args.include_role or [])
    if include_roles:
        observed_roles = {row.get(args.role_field, "") for row in manifest_rows}
        missing_roles = sorted(include_roles - observed_roles)
        if missing_roles:
            raise ValueError(f"requested comparison roles absent: {missing_roles}")
        manifest_rows = [
            row for row in manifest_rows if row.get(args.role_field, "") in include_roles
        ]
    if len(manifest_rows) != args.expected_count:
        raise ValueError(
            f"expected {args.expected_count} selected manifest rows, found {len(manifest_rows)}"
        )
    run_ids = [row["run_accession"] for row in manifest_rows]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("selected external manifest has duplicate run accessions")
    for row in manifest_rows:
        if not row.get(args.patient_id_field, "").strip():
            raise ValueError(
                f"{row['run_accession']} lacks patient field {args.patient_id_field!r}"
            )

    expression_samples, all_values = _load_expression(args.expression)
    sample_index = {sample: index for index, sample in enumerate(expression_samples)}
    missing = sorted(set(run_ids) - set(sample_index))
    if missing:
        raise ValueError(f"selected runs absent from expression matrix: {missing}")
    selected_indices = [sample_index[run] for run in run_ids]
    values = {
        gene: [row[index] for index in selected_indices]
        for gene, row in all_values.items()
    }
    zscores = _z_scores(values, run_ids)
    collectri_edges = _load_edges(args.collectri)
    pkn_edges = _load_edges(args.pkn)

    conditions: list[dict[str, object]] = []
    graph_union: set[tuple[str, str, int]] = set()
    for index, row in enumerate(manifest_rows):
        scores, target_counts = _regulon_scores(
            collectri_edges, zscores, index, args.min_targets
        )
        inputs, outputs, network = _select_network(
            pkn_edges,
            scores,
            zscores,
            index,
            args.max_outputs,
            args.max_inputs,
            args.max_depth,
            args.max_edges,
        )
        input_values = {
            node: 1.0 if zscores[node][index] >= 0 else -1.0
            for node in inputs
            if node in zscores
        }
        output_values = {node: 1.0 if value >= 0 else -1.0 for node, value in outputs}
        present_nodes = {node for source, target, _ in network for node in (source, target)}
        input_values = {
            node: value for node, value in input_values.items() if node in present_nodes
        }
        output_values = {
            node: value for node, value in output_values.items() if node in present_nodes
        }
        usable = bool(network and input_values and output_values)
        if usable:
            graph_union.update(network)
        conditions.append(
            {
                "study_accession": args.study,
                "run_accession": row["run_accession"],
                "canonical_ocm_id": row.get("canonical_ocm_id", ""),
                "patient_id": row[args.patient_id_field],
                "comparison_role": row.get(args.role_field, ""),
                "site": row.get("normalized_site", row.get("site_category", "")),
                "preprocessing_status": "included" if usable else "blocked",
                "blocked_reason": None
                if usable
                else "no_bounded_signed_path_between_expression_priors_and_TFs",
                "inputs": input_values,
                "outputs": output_values,
                "candidate_edge_count": len(network),
                "top_tf_scores": [
                    {
                        "tf": tf,
                        "score": score,
                        "target_count": target_counts.get(tf, 0),
                    }
                    for tf, score in outputs
                ],
            }
        )
    included = sum(row["preprocessing_status"] == "included" for row in conditions)
    if included < 2 or not graph_union:
        raise ValueError("fewer than two usable external conditions or empty graph union")
    result: dict[str, object] = {
        "schema_version": "regulatory_multisample_input.v1",
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
        "response_blind": True,
        "external_validation": True,
        "feature_standardization": (
            "log1p(non-negative normalized expression), then gene z-score within the "
            "selected external dataset only"
        ),
        "comparison_contract": (
            "external labels select predeclared groups only; they do not select Taylor features"
        ),
        "input_counts": {
            "selected_conditions": len(conditions),
            "included_conditions": included,
            "blocked_conditions": len(conditions) - included,
            "unique_patients": len({row["patient_id"] for row in conditions}),
            "comparison_roles": dict(
                sorted(Counter(str(row["comparison_role"]) for row in conditions).items())
            ),
            "expression_columns_total": len(expression_samples),
            "expression_genes": len(values),
            "collectri_edges": len(collectri_edges),
            "pkn_edges": len(pkn_edges),
            "frozen_graph_union_edges": len(graph_union),
        },
        "parameters": {
            "min_targets": args.min_targets,
            "max_inputs": args.max_inputs,
            "max_outputs": args.max_outputs,
            "max_depth": args.max_depth,
            "max_edges_per_condition": args.max_edges,
            "include_roles": sorted(include_roles),
            "patient_id_field": args.patient_id_field,
        },
        "sources": {
            "manifest": {"path": str(args.manifest), "sha256": _sha(args.manifest)},
            "expression": {"path": str(args.expression), "sha256": _sha(args.expression)},
            "collectri": {"path": str(args.collectri), "sha256": _sha(args.collectri)},
            "pkn": {"path": str(args.pkn), "sha256": _sha(args.pkn)},
        },
        "graph": [
            {"source": source, "target": target, "sign": sign}
            for source, target, sign in sorted(graph_union)
        ],
        "conditions": conditions,
        "claim_limit": (
            "External response-blind model input; patient-level transportability must be "
            "tested only against the already frozen Taylor signature."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=args.output.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(args.output)
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--expression", type=Path, required=True)
    result.add_argument("--manifest", type=Path, required=True)
    result.add_argument("--collectri", type=Path, required=True)
    result.add_argument("--pkn", type=Path, required=True)
    result.add_argument("--study", required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--expected-count", type=int, required=True)
    result.add_argument("--include-role", action="append")
    result.add_argument("--role-field", default="comparison_role")
    result.add_argument("--patient-id-field", default="patient_id")
    result.add_argument("--min-targets", type=int, default=5)
    result.add_argument("--max-inputs", type=int, default=3)
    result.add_argument("--max-outputs", type=int, default=6)
    result.add_argument("--max-depth", type=int, default=3)
    result.add_argument("--max-edges", type=int, default=3000)
    return result


if __name__ == "__main__":
    arguments = parser().parse_args()
    print(json.dumps(build_bundle(arguments), indent=2, sort_keys=True))
