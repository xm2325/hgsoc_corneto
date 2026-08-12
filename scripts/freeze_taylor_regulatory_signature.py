#!/usr/bin/env python3
"""Freeze a patient/cohort-recurrent Taylor regulatory edge signature.

The rule is fixed before full external inference: at nominal lambda 0.001, an
edge must occur in at least 10% of Taylor patients and in at least half of the
four cohort-specific joint fits.  Patient-balanced and richer-PKN receipts are
used only to annotate sensitivity, not to tune or exclude the primary edges.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

EXPECTED_STUDIES = {
    "E-MTAB-7223": 9,
    "E-MTAB-10801": 13,
    "E-MTAB-11000": 11,
    "E-MTAB-14568": 27,
}


class SignatureFreezeError(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_receipt(path: Path, *, expected_conditions: int) -> dict[str, Any]:
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SignatureFreezeError(f"cannot read receipt {path}: {error}") from error
    method = receipt.get("method", {})
    scope = receipt.get("scope_counts", {})
    conditions = receipt.get("conditions")
    if (
        receipt.get("status") != "completed"
        or receipt.get("response_blind") is not True
        or method.get("name") != "CarnivalFlow"
        or method.get("single_joint_problem") is not True
        or method.get("lambda_scaling") != "mean_fit"
        or float(method.get("lambda_nominal", -1)) != 0.001
        or scope.get("included_conditions") != expected_conditions
        or scope.get("preprocessing_blocked") != 0
        or not isinstance(conditions, list)
        or len(conditions) != expected_conditions
    ):
        raise SignatureFreezeError(f"receipt violates the frozen joint-lambda contract: {path}")
    for condition in conditions:
        if condition.get("status") not in {"optimal", "optimal_inaccurate"}:
            raise SignatureFreezeError(f"non-optimal Taylor condition in {path}")
        if not condition.get("run_accession") or not condition.get("patient_id"):
            raise SignatureFreezeError(f"condition lacks run/patient identity in {path}")
    return receipt


def edge_set(condition: dict[str, Any]) -> set[tuple[str, str, int]]:
    result: set[tuple[str, str, int]] = set()
    for edge in condition.get("selected_edges", []):
        source = str(edge.get("source", "")).strip()
        target = str(edge.get("target", "")).strip()
        try:
            sign = int(edge.get("sign"))
        except (TypeError, ValueError) as error:
            raise SignatureFreezeError("edge sign must be -1 or 1") from error
        if not source or not target or sign not in {-1, 1}:
            raise SignatureFreezeError("invalid selected edge")
        result.add((source, target, sign))
    return result


def atomic_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise SignatureFreezeError(f"refusing to overwrite {path}")
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise SignatureFreezeError(f"refusing to overwrite {path}")
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def freeze(args: argparse.Namespace) -> dict[str, Any]:
    pooled = read_receipt(args.pooled, expected_conditions=60)
    cohort_paths: dict[str, Path] = {}
    for item in args.cohort:
        study, separator, raw_path = item.partition("=")
        if not separator or study in cohort_paths:
            raise SignatureFreezeError(f"invalid or duplicate --cohort {item!r}")
        cohort_paths[study] = Path(raw_path)
    if set(cohort_paths) != set(EXPECTED_STUDIES):
        raise SignatureFreezeError("the four frozen Taylor cohorts are required exactly once")

    pooled_bundle = pooled.get("bundle", {}).get("sha256")
    cohort_receipts: dict[str, dict[str, Any]] = {}
    for study, expected in EXPECTED_STUDIES.items():
        receipt = read_receipt(cohort_paths[study], expected_conditions=expected)
        if receipt.get("analysis_mode") != "cohort" or receipt.get("study_accession") != study:
            raise SignatureFreezeError(f"cohort receipt mislabeled for {study}")
        if receipt.get("bundle", {}).get("sha256") != pooled_bundle:
            raise SignatureFreezeError("pooled/cohort bundle SHA-256 mismatch")
        cohort_receipts[study] = receipt

    conditions = pooled["conditions"]
    if len({row["run_accession"] for row in conditions}) != 60:
        raise SignatureFreezeError("pooled receipt does not contain 60 unique runs")
    study_counts = Counter(row.get("study_accession") for row in conditions)
    if dict(study_counts) != EXPECTED_STUDIES:
        raise SignatureFreezeError(f"unexpected pooled study counts: {dict(study_counts)}")
    patients: dict[str, set[tuple[str, str, int]]] = defaultdict(set)
    for condition in conditions:
        patients[condition["patient_id"]].update(edge_set(condition))
    if len(patients) != 52:
        raise SignatureFreezeError(f"expected 52 Taylor patients, observed {len(patients)}")
    patient_prevalence = Counter(edge for edges in patients.values() for edge in edges)
    cohort_recurrence: Counter[tuple[str, str, int]] = Counter()
    for receipt in cohort_receipts.values():
        union = set().union(*(edge_set(row) for row in receipt["conditions"]))
        cohort_recurrence.update(union)

    patient_balanced = read_receipt(args.patient_balanced, expected_conditions=52)
    richer = read_receipt(args.richer_pooled, expected_conditions=60)
    patient_balanced_union = set().union(
        *(edge_set(row) for row in patient_balanced["conditions"])
    )
    richer_union = set().union(*(edge_set(row) for row in richer["conditions"]))

    min_patients = math.ceil(args.min_patient_fraction * len(patients))
    min_cohorts = math.ceil(args.min_cohort_fraction * len(EXPECTED_STUDIES))
    selected = sorted(
        edge
        for edge, count in patient_prevalence.items()
        if count >= min_patients and cohort_recurrence[edge] >= min_cohorts
    )
    if not selected:
        raise SignatureFreezeError("the predeclared rule selected no edges")
    rows = []
    for source, target, sign in selected:
        sign_label = "+" if sign == 1 else "-"
        rows.append(
            {
                "feature_type": "edge",
                "feature_id": f"{source}|{target}|{sign_label}",
                "expected_direction": sign,
                "source": source,
                "target": target,
                "sign": sign,
                "taylor_patient_count": patient_prevalence[(source, target, sign)],
                "taylor_patient_fraction": format(
                    patient_prevalence[(source, target, sign)] / len(patients), ".12g"
                ),
                "cohort_recurrence": cohort_recurrence[(source, target, sign)],
                "patient_balanced_present": int((source, target, sign) in patient_balanced_union),
                "richer_pkn_present": int((source, target, sign) in richer_union),
            }
        )
    atomic_tsv(args.output, rows)
    receipt = {
        "schema_version": "taylor_regulatory_signature.v1",
        "status": "completed",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "response_blind": True,
        "signature_frozen_before_full_external_inference": True,
        "rule": {
            "lambda_nominal": 0.001,
            "min_patient_fraction": args.min_patient_fraction,
            "min_patient_count": min_patients,
            "min_cohort_fraction": args.min_cohort_fraction,
            "min_cohort_count": min_cohorts,
            "patient_balanced_and_richer_pkn": "annotation_only_not_selection",
        },
        "taylor_scope": {
            "runs": 60,
            "patients": 52,
            "studies": EXPECTED_STUDIES,
        },
        "selected_edges": len(rows),
        "signature": {"path": str(args.output), "sha256": sha256(args.output)},
        "inputs": {
            "pooled": {"path": str(args.pooled), "sha256": sha256(args.pooled)},
            "cohorts": {
                study: {"path": str(path), "sha256": sha256(path)}
                for study, path in sorted(cohort_paths.items())
            },
            "patient_balanced": {
                "path": str(args.patient_balanced),
                "sha256": sha256(args.patient_balanced),
            },
            "richer_pooled": {
                "path": str(args.richer_pooled),
                "sha256": sha256(args.richer_pooled),
            },
        },
        "claim_limit": (
            "Exploratory response-blind Taylor-derived regulatory signature. "
            "External datasets may test transportability but cannot retroactively tune this rule."
        ),
    }
    atomic_json(args.receipt, receipt)
    return receipt


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--pooled", type=Path, required=True)
    result.add_argument("--cohort", action="append", required=True)
    result.add_argument("--patient-balanced", type=Path, required=True)
    result.add_argument("--richer-pooled", type=Path, required=True)
    result.add_argument("--min-patient-fraction", type=float, default=0.10)
    result.add_argument("--min-cohort-fraction", type=float, default=0.50)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--receipt", type=Path, required=True)
    return result


if __name__ == "__main__":
    arguments = parser().parse_args()
    if not 0 < arguments.min_patient_fraction <= 1:
        raise SystemExit("--min-patient-fraction must be in (0,1]")
    if not 0 < arguments.min_cohort_fraction <= 1:
        raise SystemExit("--min-cohort-fraction must be in (0,1]")
    print(json.dumps(freeze(arguments), indent=2, sort_keys=True))
