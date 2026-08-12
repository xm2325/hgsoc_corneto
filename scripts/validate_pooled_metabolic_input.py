#!/usr/bin/env python3
"""Validate the frozen pooled-primary expression input for metabolic CORNETO."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any


EXPECTED_STUDIES = {
    "E-MTAB-7223": 9,
    "E-MTAB-10801": 13,
    "E-MTAB-11000": 11,
    "E-MTAB-14568": 27,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _true(value: str) -> bool:
    return value.strip().casefold() in {"1", "true", "yes"}


def validate(
    matrix: Path,
    manifest: Path,
    receipt_path: Path,
    human_gem: Path,
) -> dict[str, Any]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if not isinstance(receipt, dict) or receipt.get("status") != "completed":
        raise ValueError("pooled expression receipt is not completed")
    if receipt.get("value_name") != "tpm":
        raise ValueError("pooled matrix must contain raw TPM before log1p transformation")
    if receipt.get("primary_run_count") != 60 or receipt.get("primary_ocm_count") != 60:
        raise ValueError("pooled receipt does not contain exactly 60 primary OCM runs")
    if receipt.get("patient_count") != 52 or receipt.get("study_counts") != EXPECTED_STUDIES:
        raise ValueError("pooled receipt patient/study contract differs from the frozen cohort")

    outputs = receipt.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError("pooled receipt outputs are missing")
    matrix_record = outputs.get("matrix")
    manifest_record = outputs.get("sample_manifest")
    if not isinstance(matrix_record, dict) or not isinstance(manifest_record, dict):
        raise ValueError("pooled receipt output records are malformed")
    if _sha256(matrix) != matrix_record.get("sha256"):
        raise ValueError("pooled expression matrix SHA256 does not match its receipt")
    if _sha256(manifest) != manifest_record.get("sha256"):
        raise ValueError("pooled sample manifest SHA256 does not match its receipt")

    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 60:
        raise ValueError(f"expected 60 pooled manifest rows, found {len(rows)}")
    runs = [row["run_accession"] for row in rows]
    ocms = [row["canonical_ocm_id"] for row in rows]
    patients = [row["patient_id"] for row in rows]
    if len(set(runs)) != 60 or len(set(ocms)) != 60 or len(set(patients)) != 52:
        raise ValueError("pooled run/OCM/patient cardinality is inconsistent")
    if Counter(row["study_accession"] for row in rows) != Counter(EXPECTED_STUDIES):
        raise ValueError("pooled manifest study counts differ from the frozen cohort")
    for row in rows:
        if (
            row.get("sample_class") != "tumour"
            or row.get("histotype_group") != "HGSOC"
            or not _true(row.get("primary_cohort_eligible", ""))
        ):
            raise ValueError(f"non-primary HGSOC row in pooled manifest: {row.get('run_accession')}")

    with gzip.open(matrix, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:2] != ["gene_id", "gene_name"] or header[2:] != runs:
            raise ValueError("pooled expression columns do not exactly match manifest order")
        gene_count = sum(1 for _ in reader)
    if gene_count != receipt.get("gene_count") or gene_count != 60609:
        raise ValueError("pooled expression gene count differs from its receipt")
    if not human_gem.is_file() or human_gem.stat().st_size <= 0:
        raise ValueError("Human-GEM input is absent or empty")

    return {
        "status": "valid",
        "response_blind": True,
        "sample_count": 60,
        "ocm_count": 60,
        "patient_count": 52,
        "study_counts": EXPECTED_STUDIES,
        "gene_count": gene_count,
        "inputs": {
            "matrix": {"path": str(matrix), "sha256": _sha256(matrix)},
            "sample_manifest": {"path": str(manifest), "sha256": _sha256(manifest)},
            "pooled_expression_receipt": {
                "path": str(receipt_path),
                "sha256": _sha256(receipt_path),
            },
            "human_gem": {"path": str(human_gem), "sha256": _sha256(human_gem)},
        },
        "next_step": "Run a four-condition joint-only smoke before the 60-condition solve.",
        "claim_limit": "Input/code gate only; no metabolic solution or phenotype claim.",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--pooled-receipt", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = validate(
            args.matrix,
            args.sample_manifest,
            args.pooled_receipt,
            args.human_gem,
        )
    except (OSError, ValueError, json.JSONDecodeError, KeyError) as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"status": "valid", "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
