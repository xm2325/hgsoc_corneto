#!/usr/bin/env python3
"""Convert a formal external joint receipt into patient-level frozen-edge evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def edge_tuple(row: dict[str, Any]) -> tuple[str, str, int]:
    source, target = str(row.get("source", "")).strip(), str(row.get("target", "")).strip()
    sign = int(row.get("sign", 0))
    if not source or not target or sign not in {-1, 1}:
        raise ValueError("invalid edge record")
    return source, target, sign


def read_signature(path: Path) -> dict[tuple[str, str, int], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {edge_tuple(row): row for row in rows}
    if not rows or len(result) != len(rows):
        raise ValueError("signature is empty or contains duplicate/invalid edges")
    return result


def read_manifest(path: Path, study: str) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = {
        row["run_accession"]: row
        for row in rows
        if row.get("study_accession", study) == study
    }
    if not selected or len(selected) != sum(
        row.get("study_accession", study) == study for row in rows
    ):
        raise ValueError("manifest is empty or has duplicate runs")
    return selected


def atomic_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def build(args: argparse.Namespace) -> dict[str, Any]:
    solution, bundle = read_json(args.solution), read_json(args.bundle)
    signature = read_signature(args.signature)
    manifest = read_manifest(args.manifest, args.study)
    if (
        solution.get("status") != "completed"
        or solution.get("response_blind") is not True
        or solution.get("method", {}).get("lambda_nominal") != 0.001
        or solution.get("method", {}).get("lambda_scaling") != "mean_fit"
        or solution.get("solver", {}).get("selected") != "GUROBI"
        or solution.get("solver", {}).get("has_incumbent") is not True
    ):
        raise ValueError("formal external solution violates the frozen solver contract")
    if solution.get("bundle", {}).get("sha256") != sha256(args.bundle):
        raise ValueError("solution/bundle SHA-256 mismatch")
    required = bundle.get("sources", {}).get("required_signature", {})
    if required.get("sha256") != sha256(args.signature):
        raise ValueError("bundle was not constructed with this frozen signature")
    candidate_edges = {edge_tuple(row) for row in bundle.get("graph", [])}
    if not set(signature).issubset(candidate_edges):
        raise ValueError("not every frozen signature edge is in the candidate graph")
    conditions = solution.get("conditions", [])
    if not conditions or any(
        row.get("status") not in {"optimal", "optimal_inaccurate"} for row in conditions
    ):
        raise ValueError("external solution has missing or non-optimal conditions")
    by_run = {row["run_accession"]: row for row in conditions}
    if len(by_run) != len(conditions):
        raise ValueError("solution has duplicate run accessions")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {"groups": {}}
    seen_labels: set[str] = set()
    for spec in args.group:
        label, separator, rule = spec.partition("=")
        field, second_separator, value = rule.partition(":")
        if (
            not separator
            or not second_separator
            or not re.fullmatch(r"[A-Za-z0-9_]+", label)
            or label in seen_labels
        ):
            raise ValueError(f"invalid group specification {spec!r}")
        seen_labels.add(label)
        selected_runs = [
            run for run, row in manifest.items() if row.get(field, "") == value and run in by_run
        ]
        if not selected_runs:
            raise ValueError(f"group {label} selected no solved conditions")
        patient_conditions: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for run in selected_runs:
            patient = manifest[run].get(args.patient_id_field, "").strip()
            if not patient:
                raise ValueError(f"{run} lacks patient identifier")
            patient_conditions[patient].append(by_run[run])
        if len(patient_conditions) < 2:
            raise ValueError(f"group {label} has fewer than two patients")
        evidence_rows: list[dict[str, Any]] = []
        for patient in sorted(patient_conditions):
            rows = patient_conditions[patient]
            selected_by_condition = [
                {edge_tuple(edge) for edge in row.get("selected_edges", [])} for row in rows
            ]
            for edge, signature_row in sorted(signature.items()):
                selected_count = sum(edge in selected for selected in selected_by_condition)
                prevalence = selected_count / len(rows)
                selected = prevalence >= args.within_patient_threshold
                evidence_rows.append(
                    {
                        "patient_id": patient,
                        "feature_type": "edge",
                        "feature_id": signature_row["feature_id"],
                        "evaluable": 1,
                        "selected": int(selected),
                        "direction": edge[2] if selected else 0,
                        "condition_count": len(rows),
                        "selected_condition_count": selected_count,
                        "condition_prevalence": format(prevalence, ".12g"),
                    }
                )
        evidence_path = args.output_dir / f"{label}.evidence.tsv"
        contract_path = args.output_dir / f"{label}.contract.json"
        if evidence_path.exists() or contract_path.exists():
            raise ValueError(f"refusing to overwrite group outputs for {label}")
        atomic_tsv(evidence_path, evidence_rows)
        contract = {
            "schema_version": "external_corneto_group.v1",
            "status": "completed",
            "created_at_utc": datetime.now(UTC).isoformat(),
            "group_label": label,
            "source_accession": args.study,
            "evidence_sha256": sha256(evidence_path),
            "signature_sha256": sha256(args.signature),
            "patient_count": len(patient_conditions),
            "analysis_unit": args.analysis_unit,
            "normalization": {
                "performed_within_dataset": True,
                "pooled_raw_expression": False,
                "input_scale": "network_selection",
                "receipt_path": str(args.normalization_receipt),
                "receipt_sha256": sha256(args.normalization_receipt),
            },
            "independence": {
                "cells_as_replicates": False,
                "patient_id_column": "patient_id",
                "within_patient_condition_aggregation": (
                    f"selected_prevalence >= {args.within_patient_threshold}"
                ),
            },
            "inference": {
                "signature_frozen_before_external_scoring": True,
                "feature_selection_using_external_labels": False,
                "all_frozen_edges_present_in_candidate_graph": True,
            },
            "group_rule": {"field": field, "value": value},
            "solution": {"path": str(args.solution), "sha256": sha256(args.solution)},
            "bundle": {"path": str(args.bundle), "sha256": sha256(args.bundle)},
            "claim_limit": (
                "Patient-level model-selection evidence; not clinical response, causality, "
                "or measured signalling activity."
            ),
        }
        atomic_json(contract_path, contract)
        result["groups"][label] = {
            "patients": len(patient_conditions),
            "conditions": len(selected_runs),
            "evidence": str(evidence_path),
            "contract": str(contract_path),
        }
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--solution", type=Path, required=True)
    result.add_argument("--bundle", type=Path, required=True)
    result.add_argument("--signature", type=Path, required=True)
    result.add_argument("--manifest", type=Path, required=True)
    result.add_argument("--normalization-receipt", type=Path, required=True)
    result.add_argument("--study", required=True)
    result.add_argument("--patient-id-field", required=True)
    result.add_argument(
        "--analysis-unit",
        choices=("patient_model", "patient_pseudobulk", "patient_tissue"),
        required=True,
    )
    result.add_argument("--group", action="append", required=True)
    result.add_argument("--within-patient-threshold", type=float, default=0.5)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


if __name__ == "__main__":
    arguments = parser().parse_args()
    if not 0 < arguments.within_patient_threshold <= 1:
        raise SystemExit("--within-patient-threshold must be in (0,1]")
    print(json.dumps(build(arguments), indent=2, sort_keys=True))
