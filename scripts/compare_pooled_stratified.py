#!/usr/bin/env python3
"""Compare pooled and cohort-stratified regulatory/NMF results.

The script is deliberately solver-free.  It joins samples through the frozen
manifest, compares regulatory edge sets at each lambda, and aligns cohort-local
NMF states to pooled states by maximum sample overlap.  Local NMF state labels
are never assumed to have the same meaning across studies.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable


Edge = tuple[str, str, int]


class ComparisonError(ValueError):
    """Raised when an input violates the comparison contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _lambda_label(value: str | float) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ComparisonError(f"invalid lambda value {value!r}") from error
    if not math.isfinite(number) or number < 0:
        raise ComparisonError(f"invalid lambda value {value!r}")
    return f"{number:.12g}"


def _jaccard(left: set[Any], right: set[Any]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def _read_manifest(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"run_accession", "study_accession", "canonical_ocm_id", "patient_id"}
    if not rows or not required.issubset(rows[0]):
        raise ComparisonError("manifest is empty or missing required columns")
    by_run: dict[str, dict[str, str]] = {}
    for line, row in enumerate(rows, start=2):
        run = row["run_accession"].strip()
        if not run or run in by_run:
            raise ComparisonError(f"manifest line {line} has invalid/duplicate run_accession")
        by_run[run] = row
    primary = [row for row in rows if _truthy(row.get("primary_cohort_eligible", ""))]
    summary = {
        "all_run_count": len(rows),
        "primary_run_count": len(primary),
        "primary_ocm_count": len({row["canonical_ocm_id"] for row in primary}),
        "primary_patient_count": len({row["patient_id"] for row in primary}),
        "primary_runs_by_study": dict(sorted(Counter(row["study_accession"] for row in primary).items())),
    }
    return by_run, summary


def _edge(edge: Any, where: str) -> Edge:
    if not isinstance(edge, dict):
        raise ComparisonError(f"{where} is not an object")
    source, target, sign = edge.get("source"), edge.get("target"), edge.get("sign")
    if not isinstance(source, str) or not source or not isinstance(target, str) or not target:
        raise ComparisonError(f"{where} has an invalid source or target")
    if sign not in (-1, 1):
        raise ComparisonError(f"{where} sign must be -1 or 1")
    return source, target, int(sign)


def _declares_joint(root: dict[str, Any]) -> bool:
    method = root.get("method") if isinstance(root.get("method"), dict) else {}
    candidates = (
        root.get("joint_multi_sample"),
        root.get("joint_inference"),
        method.get("joint_multi_sample"),
        method.get("joint_inference"),
    )
    return any(value is True for value in candidates)


def _read_regulatory(
    path: Path,
    manifest: dict[str, dict[str, str]],
    expected_lambda: str,
) -> dict[str, Any]:
    try:
        root = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ComparisonError(f"cannot read regulatory receipt {path}: {error}") from error
    if not isinstance(root, dict) or root.get("status") != "completed":
        raise ComparisonError(f"regulatory receipt {path} must have status='completed'")
    samples = root.get("samples")
    if not isinstance(samples, list) or not samples:
        raise ComparisonError(f"regulatory receipt {path} has no samples")
    method = root.get("method") if isinstance(root.get("method"), dict) else {}
    observed_lambda = method.get(
        "lambda_reg_reported", method.get("lambda_reg", root.get("lambda_reg"))
    )
    if observed_lambda is not None and _lambda_label(observed_lambda) != expected_lambda:
        raise ComparisonError(
            f"regulatory receipt {path} lambda={observed_lambda!r}, expected {expected_lambda}"
        )
    by_run: dict[str, dict[str, Any]] = {}
    for index, sample in enumerate(samples):
        if not isinstance(sample, dict):
            raise ComparisonError(f"{path}: samples[{index}] is not an object")
        run = sample.get("run_accession")
        if not isinstance(run, str) or not run or run in by_run:
            raise ComparisonError(f"{path}: samples[{index}] has invalid/duplicate run_accession")
        if run not in manifest:
            raise ComparisonError(f"{path}: run {run} is absent from manifest")
        raw_edges = sample.get("selected_edges")
        if not isinstance(raw_edges, list):
            raise ComparisonError(f"{path}: sample {run} has no selected_edges array")
        edges = {_edge(value, f"{path}: sample {run} edge") for value in raw_edges}
        if len(edges) != len(raw_edges):
            raise ComparisonError(f"{path}: sample {run} has duplicate selected edges")
        by_run[run] = {
            "study": manifest[run]["study_accession"],
            "ocm": manifest[run]["canonical_ocm_id"],
            "patient": manifest[run]["patient_id"],
            "status": str(sample.get("status", "unknown")),
            "edges": edges,
        }
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "schema_version": root.get("schema_version"),
        "declares_joint_multi_sample": _declares_joint(root),
        "samples": by_run,
    }


def _parse_pooled_regulatory(
    specs: Iterable[str], manifest: dict[str, dict[str, str]]
) -> dict[str, dict[str, Any]]:
    result = {}
    for spec in specs:
        if "=" not in spec:
            raise ComparisonError("--pooled-regulatory must be LAMBDA=JSON")
        raw_lambda, raw_path = spec.split("=", 1)
        label = _lambda_label(raw_lambda)
        if label in result:
            raise ComparisonError(f"duplicate pooled lambda {label}")
        result[label] = _read_regulatory(Path(raw_path), manifest, label)
    return result


def _parse_cohort_regulatory(
    specs: Iterable[str], manifest: dict[str, dict[str, str]]
) -> dict[str, dict[str, dict[str, Any]]]:
    result: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for spec in specs:
        fields = spec.split("|", 2)
        if len(fields) != 3:
            raise ComparisonError("--cohort-regulatory must be STUDY|LAMBDA|JSON")
        study, raw_lambda, raw_path = fields
        label = _lambda_label(raw_lambda)
        if label in result[study]:
            raise ComparisonError(f"duplicate cohort receipt {study}/{label}")
        receipt = _read_regulatory(Path(raw_path), manifest, label)
        wrong = sorted(run for run, row in receipt["samples"].items() if row["study"] != study)
        if wrong:
            raise ComparisonError(f"{study}/{label} contains runs from other studies: {wrong[:3]}")
        result[study][label] = receipt
    return dict(result)


def _regulatory_comparison(
    pooled: dict[str, dict[str, Any]],
    cohorts: dict[str, dict[str, dict[str, Any]]],
    require_joint: bool,
) -> dict[str, Any]:
    if not pooled or not cohorts:
        return {}
    if require_joint:
        invalid = sorted(label for label, receipt in pooled.items() if not receipt["declares_joint_multi_sample"])
        if invalid:
            raise ComparisonError(
                "pooled regulatory receipts do not declare joint_multi_sample=true for lambda: "
                + ", ".join(invalid)
            )
    pooled_lambdas = set(pooled)
    common_lambdas = pooled_lambdas.intersection(*(set(values) for values in cohorts.values()))
    if not common_lambdas:
        raise ComparisonError("no lambda value is shared by pooled and all cohort receipts")
    lambda_rows: list[dict[str, Any]] = []
    study_rows: list[dict[str, Any]] = []
    for label in sorted(common_lambdas, key=float):
        pooled_samples = pooled[label]["samples"]
        cohort_samples: dict[str, dict[str, Any]] = {}
        for study in sorted(cohorts):
            for run, row in cohorts[study][label]["samples"].items():
                if run in cohort_samples:
                    raise ComparisonError(f"run {run} is duplicated across cohort receipts at lambda {label}")
                cohort_samples[run] = row
        matched = sorted(set(pooled_samples) & set(cohort_samples))
        pooled_union = set().union(*(row["edges"] for row in pooled_samples.values()))
        cohort_union = set().union(*(row["edges"] for row in cohort_samples.values()))
        sample_jaccards = [
            _jaccard(pooled_samples[run]["edges"], cohort_samples[run]["edges"])
            for run in matched
        ]
        lambda_rows.append(
            {
                "lambda": float(label),
                "pooled_sample_count": len(pooled_samples),
                "cohort_sample_count": len(cohort_samples),
                "matched_sample_count": len(matched),
                "pooled_only_count": len(set(pooled_samples) - set(cohort_samples)),
                "cohort_only_count": len(set(cohort_samples) - set(pooled_samples)),
                "pooled_edge_union_size": len(pooled_union),
                "cohort_merged_edge_union_size": len(cohort_union),
                "edge_union_intersection_size": len(pooled_union & cohort_union),
                "edge_union_jaccard": _jaccard(pooled_union, cohort_union),
                "mean_matched_sample_edge_jaccard": (
                    sum(sample_jaccards) / len(sample_jaccards) if sample_jaccards else None
                ),
                "pooled_declares_joint_multi_sample": pooled[label]["declares_joint_multi_sample"],
            }
        )
        for study in sorted(cohorts):
            pooled_study = {run: row for run, row in pooled_samples.items() if row["study"] == study}
            local = cohorts[study][label]["samples"]
            matched_study = sorted(set(pooled_study) & set(local))
            pooled_study_union = set().union(*(row["edges"] for row in pooled_study.values()))
            local_union = set().union(*(row["edges"] for row in local.values()))
            per_sample = [
                _jaccard(pooled_study[run]["edges"], local[run]["edges"])
                for run in matched_study
            ]
            study_rows.append(
                {
                    "study_accession": study,
                    "lambda": float(label),
                    "pooled_sample_count": len(pooled_study),
                    "cohort_sample_count": len(local),
                    "matched_sample_count": len(matched_study),
                    "pooled_edge_union_size": len(pooled_study_union),
                    "cohort_edge_union_size": len(local_union),
                    "edge_union_intersection_size": len(pooled_study_union & local_union),
                    "edge_union_jaccard": _jaccard(pooled_study_union, local_union),
                    "mean_matched_sample_edge_jaccard": (
                        sum(per_sample) / len(per_sample) if per_sample else None
                    ),
                }
            )
    return {
        "common_lambda_values": [float(value) for value in sorted(common_lambdas, key=float)],
        "pooled_receipts": {
            label: {key: value for key, value in receipt.items() if key != "samples"}
            for label, receipt in sorted(pooled.items(), key=lambda item: float(item[0]))
        },
        "cohort_receipts": {
            study: {
                label: {key: value for key, value in receipt.items() if key != "samples"}
                for label, receipt in sorted(values.items(), key=lambda item: float(item[0]))
            }
            for study, values in sorted(cohorts.items())
        },
        "lambda_comparisons": lambda_rows,
        "study_comparisons": study_rows,
        "interpretation_gate": (
            "pooled-vs-stratified lambda interpretation is valid only when each pooled receipt "
            "declares a single joint multi-sample optimization"
        ),
    }


def _read_assignments(
    path: Path,
    manifest: dict[str, dict[str, str]],
    state_column: str | None = None,
) -> dict[str, Any]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "run_accession" not in rows[0]:
        raise ComparisonError(f"NMF assignments {path} are empty or lack run_accession")
    candidates = [state_column] if state_column else [
        "pooled_state", "independent_state", "technical_state", "nmf_state", "state", "cluster"
    ]
    selected = next((column for column in candidates if column and column in rows[0]), None)
    if selected is None:
        raise ComparisonError(f"NMF assignments {path} have no recognized state column")
    by_run: dict[str, str] = {}
    for line, row in enumerate(rows, start=2):
        run, state = row["run_accession"].strip(), row[selected].strip()
        if not run or not state or run in by_run:
            raise ComparisonError(f"NMF assignments {path} line {line} has invalid/duplicate values")
        if run not in manifest:
            raise ComparisonError(f"NMF assignments {path}: run {run} is absent from manifest")
        by_run[run] = state
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "state_column": selected,
        "assignments": by_run,
    }


def _best_state_mapping(
    cohort: dict[str, str], pooled: dict[str, str], matched: list[str]
) -> dict[str, str | None]:
    local_states = sorted({cohort[run] for run in matched})
    pooled_states = sorted({pooled[run] for run in matched})
    size = max(len(local_states), len(pooled_states))
    if size > 16:
        raise ComparisonError("more than 16 NMF states; exact state alignment is intentionally bounded")
    weights = [
        [sum(cohort[run] == left and pooled[run] == right for run in matched) for right in pooled_states]
        + [0] * (size - len(pooled_states))
        for left in local_states
    ]
    weights.extend([[0] * size for _ in range(size - len(local_states))])
    # Maximum-weight square assignment using a bit-mask dynamic program.  NMF
    # ranks are small, so this avoids a heavy SciPy dependency in the audit job.
    scores: dict[int, tuple[int, tuple[int, ...]]] = {0: (0, ())}
    for row in range(size):
        nxt: dict[int, tuple[int, tuple[int, ...]]] = {}
        for mask, (score, chosen) in scores.items():
            for column in range(size):
                if mask & (1 << column):
                    continue
                candidate = (score + weights[row][column], chosen + (column,))
                new_mask = mask | (1 << column)
                if new_mask not in nxt or candidate[0] > nxt[new_mask][0]:
                    nxt[new_mask] = candidate
        scores = nxt
    chosen = scores[(1 << size) - 1][1]
    return {
        state: pooled_states[chosen[index]] if chosen[index] < len(pooled_states) else None
        for index, state in enumerate(local_states)
    }


def _comb2(value: int) -> int:
    return value * (value - 1) // 2


def _adjusted_rand(left: dict[str, str], right: dict[str, str], matched: list[str]) -> float:
    if len(matched) < 2:
        return 1.0
    left_counts = Counter(left[run] for run in matched)
    right_counts = Counter(right[run] for run in matched)
    cells = Counter((left[run], right[run]) for run in matched)
    sum_cells = sum(_comb2(value) for value in cells.values())
    sum_left = sum(_comb2(value) for value in left_counts.values())
    sum_right = sum(_comb2(value) for value in right_counts.values())
    pairs = _comb2(len(matched))
    expected = sum_left * sum_right / pairs
    maximum = (sum_left + sum_right) / 2
    denominator = maximum - expected
    return 1.0 if denominator == 0 and sum_cells == maximum else (sum_cells - expected) / denominator


def _nmf_comparison(
    pooled: dict[str, Any] | None,
    cohorts: dict[str, dict[str, Any]],
    manifest: dict[str, dict[str, str]],
) -> dict[str, Any]:
    if pooled is None or not cohorts:
        return {}
    pooled_assignments = pooled["assignments"]
    summaries: list[dict[str, Any]] = []
    alignments: list[dict[str, Any]] = []
    for study, receipt in sorted(cohorts.items()):
        local = receipt["assignments"]
        wrong = sorted(run for run in local if manifest[run]["study_accession"] != study)
        if wrong:
            raise ComparisonError(f"NMF cohort {study} contains runs from other studies: {wrong[:3]}")
        pooled_study = {
            run: state for run, state in pooled_assignments.items()
            if manifest[run]["study_accession"] == study
        }
        matched = sorted(set(local) & set(pooled_study))
        if not matched:
            raise ComparisonError(f"NMF cohort {study} has no overlap with pooled assignments")
        mapping = _best_state_mapping(local, pooled_study, matched)
        for local_state in sorted(set(local[run] for run in matched)):
            pooled_state = mapping[local_state]
            local_runs = {run for run in matched if local[run] == local_state}
            pooled_runs = {
                run for run in matched if pooled_state is not None and pooled_study[run] == pooled_state
            }
            alignments.append(
                {
                    "study_accession": study,
                    "cohort_state": local_state,
                    "aligned_pooled_state": pooled_state,
                    "cohort_state_sample_count": len(local_runs),
                    "pooled_state_sample_count_within_study": len(pooled_runs),
                    "intersection_count": len(local_runs & pooled_runs),
                    "union_count": len(local_runs | pooled_runs),
                    "state_jaccard": _jaccard(local_runs, pooled_runs),
                    "precision_to_pooled_state": (
                        len(local_runs & pooled_runs) / len(local_runs) if local_runs else None
                    ),
                    "recall_of_pooled_state": (
                        len(local_runs & pooled_runs) / len(pooled_runs) if pooled_runs else None
                    ),
                }
            )
        summaries.append(
            {
                "study_accession": study,
                "pooled_sample_count_within_study": len(pooled_study),
                "cohort_sample_count": len(local),
                "matched_sample_count": len(matched),
                "pooled_only_count": len(set(pooled_study) - set(local)),
                "cohort_only_count": len(set(local) - set(pooled_study)),
                "pooled_state_count_within_study": len({pooled_study[run] for run in matched}),
                "cohort_state_count": len({local[run] for run in matched}),
                "adjusted_rand_index": _adjusted_rand(local, pooled_study, matched),
            }
        )
    return {
        "pooled_input": {key: value for key, value in pooled.items() if key != "assignments"},
        "cohort_inputs": {
            study: {key: value for key, value in receipt.items() if key != "assignments"}
            for study, receipt in sorted(cohorts.items())
        },
        "pooled_sample_count": len(pooled_assignments),
        "pooled_state_counts": dict(sorted(Counter(pooled_assignments.values()).items())),
        "study_summaries": summaries,
        "state_alignment": alignments,
        "alignment_method": "maximum sample-overlap bipartite assignment; labels are cohort-local",
    }


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def compare(args: argparse.Namespace) -> dict[str, Any]:
    manifest, manifest_summary = _read_manifest(args.manifest)
    pooled_regulatory = _parse_pooled_regulatory(args.pooled_regulatory, manifest)
    cohort_regulatory = _parse_cohort_regulatory(args.cohort_regulatory, manifest)
    pooled_nmf = (
        _read_assignments(args.pooled_nmf, manifest, args.pooled_state_column)
        if args.pooled_nmf else None
    )
    cohort_nmf: dict[str, dict[str, Any]] = {}
    for spec in args.cohort_nmf:
        fields = spec.split("|", 1)
        if len(fields) != 2:
            raise ComparisonError("--cohort-nmf must be STUDY|TSV")
        study, raw_path = fields
        if study in cohort_nmf:
            raise ComparisonError(f"duplicate cohort NMF assignment for {study}")
        cohort_nmf[study] = _read_assignments(Path(raw_path), manifest, args.cohort_state_column)
    regulatory = _regulatory_comparison(
        pooled_regulatory, cohort_regulatory, args.require_joint_pooled
    )
    nmf = _nmf_comparison(pooled_nmf, cohort_nmf, manifest)
    if not regulatory and not nmf:
        raise ComparisonError("provide both pooled and cohort inputs for regulatory and/or NMF")
    return {
        "schema_version": "pooled_stratified_comparison.v1",
        "status": "completed",
        "response_blind": True,
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "manifest": {
            "path": str(args.manifest),
            "sha256": _sha256(args.manifest),
            **manifest_summary,
        },
        "regulatory": regulatory,
        "nmf": nmf,
        "claim_limit": (
            "response-blind technical stability and cross-study comparison; NMF state labels are "
            "not Barnes subtypes and no drug-response or causal claim is made"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--pooled-regulatory", action="append", default=[], metavar="LAMBDA=JSON")
    parser.add_argument(
        "--cohort-regulatory", action="append", default=[], metavar="STUDY|LAMBDA|JSON"
    )
    parser.add_argument("--require-joint-pooled", action="store_true")
    parser.add_argument("--pooled-nmf", type=Path)
    parser.add_argument("--cohort-nmf", action="append", default=[], metavar="STUDY|TSV")
    parser.add_argument("--pooled-state-column")
    parser.add_argument("--cohort-state-column")
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        result = compare(args)
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output_json.with_suffix(args.output_json.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output_json)
        output_dir = args.output_dir or args.output_json.parent / "pooled_stratified_tables"
        regulatory = result["regulatory"]
        nmf = result["nmf"]
        _write_tsv(output_dir / "regulatory_lambda_comparison.tsv", regulatory.get("lambda_comparisons", []))
        _write_tsv(output_dir / "regulatory_study_comparison.tsv", regulatory.get("study_comparisons", []))
        _write_tsv(output_dir / "nmf_study_summary.tsv", nmf.get("study_summaries", []))
        _write_tsv(output_dir / "nmf_state_alignment.tsv", nmf.get("state_alignment", []))
    except (ComparisonError, OSError, json.JSONDecodeError) as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": "completed",
                "regulatory_lambdas": len(result["regulatory"].get("lambda_comparisons", [])),
                "nmf_studies": len(result["nmf"].get("study_summaries", [])),
                "output": str(args.output_json),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
