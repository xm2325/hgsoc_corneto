#!/usr/bin/env python3
"""Compare pooled-60 and one-OCM-per-patient NMF assignments."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.optimize import linear_sum_assignment


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {}
    for row in rows:
        run = row.get("run_accession", "")
        state = row.get("state", row.get("technical_state", row.get("independent_state", "")))
        if not run or not state or run in result:
            raise ValueError(f"invalid or duplicate assignment in {path}")
        result[run] = state
    if not result:
        raise ValueError(f"empty assignments: {path}")
    return result


def compare(pooled: Path, balanced: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise ValueError(f"refusing to overwrite {output}")
    left, right = _read(pooled), _read(balanced)
    common = sorted(set(left) & set(right))
    if len(common) != 52:
        raise ValueError(f"expected 52 common patient-balanced runs, found {len(common)}")
    left_states = sorted(set(left[run] for run in common))
    right_states = sorted(set(right[run] for run in common))
    matrix = [[sum(left[run] == a and right[run] == b for run in common) for b in right_states] for a in left_states]
    rows, cols = linear_sum_assignment([[-value for value in row] for row in matrix])
    mapping = {right_states[col]: left_states[row] for row, col in zip(rows, cols, strict=True)}
    mapped = [mapping.get(right[run], "unmapped") for run in common]
    pooled_labels = [left[run] for run in common]
    agreement = sum(a == b for a, b in zip(pooled_labels, mapped, strict=True)) / len(common)
    result = {
        "status": "completed",
        "schema_version": "patient_balanced_nmf_comparison.v1",
        "pooled_input": {"path": str(pooled), "sha256": _sha(pooled)},
        "patient_balanced_input": {"path": str(balanced), "sha256": _sha(balanced)},
        "common_run_count": len(common),
        "pooled_state_counts": {state: pooled_labels.count(state) for state in left_states},
        "patient_balanced_state_counts": {state: mapped.count(state) for state in sorted(set(mapped))},
        "state_mapping_patient_to_pooled": mapping,
        "mapped_assignment_agreement": agreement,
        "adjusted_rand_index": adjusted_rand_score(pooled_labels, mapped),
        "normalized_mutual_information": normalized_mutual_info_score(pooled_labels, mapped),
        "confusion_matrix_pooled_rows_patient_columns": matrix,
        "claim_limit": "response-blind patient-weighting sensitivity; no phenotype or causal interpretation",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pooled", type=Path, required=True)
    parser.add_argument("--patient-balanced", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.pooled, args.patient_balanced, args.output)
    print(json.dumps({"status": result["status"], "common_run_count": result["common_run_count"], "ari": result["adjusted_rand_index"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
