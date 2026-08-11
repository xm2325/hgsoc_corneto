#!/usr/bin/env python3
"""Freeze one auditable 60-OCM input bundle for multi-sample CARNIVAL.

The legacy regulatory pilot constructs a different bounded graph and solves a
separate one-condition problem for every OCM.  That is useful as an independent
baseline, but it cannot test CORNETO's structured multi-sample penalty.  This
script derives the same response-blind features once, within each RNA study,
then freezes their all-study graph union so pooled and cohort-stratified solves
use identical per-condition evidence and the same PKN search space.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import UTC, datetime
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

EXPECTED_COUNTS = {
    "E-MTAB-7223": 9,
    "E-MTAB-10801": 13,
    "E-MTAB-11000": 11,
    "E-MTAB-14568": 27,
}


def _matrix_spec(value: str) -> tuple[str, Path]:
    study, separator, raw_path = value.partition("=")
    if not separator or not study or not raw_path:
        raise argparse.ArgumentTypeError("expression must be STUDY=PATH")
    return study, Path(raw_path)


def build_bundle(
    *,
    expressions: list[tuple[str, Path]],
    manifest: Path,
    collectri: Path,
    pkn: Path,
    output: Path,
    min_targets: int,
    max_inputs: int,
    max_outputs: int,
    max_depth: int,
    max_edges: int,
) -> dict[str, object]:
    studies = [study for study, _ in expressions]
    if len(studies) != len(set(studies)):
        raise ValueError("duplicate study expression arguments")
    if set(studies) != set(EXPECTED_COUNTS):
        raise ValueError(
            f"expected studies {sorted(EXPECTED_COUNTS)}, received {sorted(studies)}"
        )
    for path in (manifest, collectri, pkn, *(path for _, path in expressions)):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output.exists():
        raise ValueError(f"refusing to overwrite {output}")

    collectri_edges = _load_edges(collectri)
    pkn_edges = _load_edges(pkn)
    conditions: list[dict[str, object]] = []
    graph_union: set[tuple[str, str, int]] = set()
    seen_runs: set[str] = set()
    cohort_expression_counts: dict[str, int] = {}

    # Scores are standardized within study, but the resulting signed features
    # are frozen before any pooled/cohort solve.  Thus the two formulations see
    # identical evidence and cannot differ because of re-standardization.
    for study, expression_path in expressions:
        manifest_rows = _load_manifest(manifest, study, primary_only=True)
        expected = EXPECTED_COUNTS[study]
        if len(manifest_rows) != expected:
            raise ValueError(
                f"{study}: expected {expected} primary rows, found {len(manifest_rows)}"
            )
        expression_samples, values = _load_expression(expression_path)
        cohort_expression_counts[study] = len(expression_samples)
        sample_index = {sample: index for index, sample in enumerate(expression_samples)}
        missing = sorted(
            row["run_accession"]
            for row in manifest_rows
            if row["run_accession"] not in sample_index
        )
        if missing:
            raise ValueError(f"{study}: primary runs absent from expression matrix: {missing}")
        zscores = _z_scores(values, expression_samples)
        for row in manifest_rows:
            run = row["run_accession"]
            if run in seen_runs:
                raise ValueError(f"duplicate primary run across studies: {run}")
            seen_runs.add(run)
            index = sample_index[run]
            scores, target_counts = _regulon_scores(collectri_edges, zscores, index, min_targets)
            inputs, outputs, network = _select_network(
                pkn_edges,
                scores,
                zscores,
                index,
                max_outputs,
                max_inputs,
                max_depth,
                max_edges,
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
                    "study_accession": study,
                    "run_accession": run,
                    "canonical_ocm_id": row.get("canonical_ocm_id", ""),
                    "patient_id": row.get("patient_id", ""),
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

    if len(conditions) != 60 or len(seen_runs) != 60:
        raise ValueError(f"expected 60 unique primary conditions, found {len(seen_runs)}")
    if not graph_union:
        raise ValueError("all-primary graph union is empty")
    study_counts = Counter(str(row["study_accession"]) for row in conditions)
    included_counts = Counter(
        str(row["study_accession"])
        for row in conditions
        if row["preprocessing_status"] == "included"
    )
    result: dict[str, object] = {
        "schema_version": "regulatory_multisample_input.v1",
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
        "response_blind": True,
        "feature_standardization": "log1p(TPM) gene z-score within RNA study",
        "comparison_contract": (
            "pooled and cohort-stratified solvers must use this same graph and these exact "
            "per-condition signed features"
        ),
        "input_counts": {
            "primary_conditions": len(conditions),
            "included_conditions": sum(included_counts.values()),
            "blocked_conditions": len(conditions) - sum(included_counts.values()),
            "study_conditions": dict(sorted(study_counts.items())),
            "study_included_conditions": dict(sorted(included_counts.items())),
            "study_expression_columns": cohort_expression_counts,
            "collectri_edges": len(collectri_edges),
            "pkn_edges": len(pkn_edges),
            "frozen_graph_union_edges": len(graph_union),
        },
        "parameters": {
            "min_targets": min_targets,
            "max_inputs": max_inputs,
            "max_outputs": max_outputs,
            "max_depth": max_depth,
            "max_edges_per_condition": max_edges,
        },
        "sources": {
            "manifest": {"path": str(manifest), "sha256": _sha(manifest)},
            "collectri": {"path": str(collectri), "sha256": _sha(collectri)},
            "pkn": {"path": str(pkn), "sha256": _sha(pkn)},
            "expressions": [
                {"study_accession": study, "path": str(path), "sha256": _sha(path)}
                for study, path in expressions
            ],
        },
        "graph": [
            {"source": source, "target": target, "sign": sign}
            for source, target, sign in sorted(graph_union)
        ],
        "conditions": conditions,
        "claim_limit": "response-blind method input; no drug-response or causal interpretation",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", action="append", type=_matrix_spec, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--collectri", type=Path, required=True)
    parser.add_argument("--pkn", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-targets", type=int, default=5)
    parser.add_argument("--max-inputs", type=int, default=3)
    parser.add_argument("--max-outputs", type=int, default=6)
    parser.add_argument("--max-depth", type=int, default=3)
    parser.add_argument("--max-edges", type=int, default=3000)
    args = parser.parse_args()
    result = build_bundle(
        expressions=args.expression,
        manifest=args.manifest,
        collectri=args.collectri,
        pkn=args.pkn,
        output=args.output,
        min_targets=args.min_targets,
        max_inputs=args.max_inputs,
        max_outputs=args.max_outputs,
        max_depth=args.max_depth,
        max_edges=args.max_edges,
    )
    print(json.dumps({"status": result["status"], **result["input_counts"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
