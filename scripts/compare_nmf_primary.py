#!/usr/bin/env python3
"""Compare pooled-60 and cohort-stratified primary-OCM NMF at fixed rank 3."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import chi2_contingency
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _aligned_agreement(reference: list[str], candidate: list[str]) -> tuple[float, dict[str, str]]:
    ref_names = sorted(set(reference))
    cand_names = sorted(set(candidate))
    table = np.zeros((len(cand_names), len(ref_names)), dtype=int)
    for candidate_label, reference_label in zip(candidate, reference, strict=True):
        table[cand_names.index(candidate_label), ref_names.index(reference_label)] += 1
    rows, columns = linear_sum_assignment(-table)
    mapping = {cand_names[row]: ref_names[column] for row, column in zip(rows, columns, strict=True)}
    agreement = np.mean([mapping.get(value) == expected for value, expected in zip(candidate, reference, strict=True)])
    return float(agreement), mapping


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError(f"Refusing to overwrite {args.output_dir}")
    pooled_path = args.root / "pooled60" / "rank_3" / "assignments.tsv"
    pooled = _rows(pooled_path)
    if len(pooled) != 60:
        raise ValueError(f"Expected 60 pooled primary OCMs, found {len(pooled)}")
    pooled_by_run = {row["run_accession"]: row for row in pooled}
    studies = ("E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568")
    comparison: list[dict[str, object]] = []
    for study in studies:
        cohort = _rows(args.root / study / "rank_3" / "assignments.tsv")
        expected = [pooled_by_run[row["run_accession"]]["technical_state"] for row in cohort]
        observed = [row["technical_state"] for row in cohort]
        agreement, mapping = _aligned_agreement(expected, observed)
        comparison.append(
            {
                "study_accession": study,
                "sample_count": len(cohort),
                "adjusted_rand_index": format(adjusted_rand_score(expected, observed), ".12g"),
                "normalized_mutual_information": format(normalized_mutual_info_score(expected, observed), ".12g"),
                "hungarian_aligned_agreement": format(agreement, ".12g"),
                "cohort_to_pooled_state_mapping": json.dumps(mapping, sort_keys=True),
            }
        )
    study_names = sorted({row["study_accession"] for row in pooled})
    state_names = sorted({row["technical_state"] for row in pooled})
    table = np.asarray(
        [[sum(row["study_accession"] == study and row["technical_state"] == state for row in pooled) for state in state_names] for study in study_names],
        dtype=int,
    )
    chi_square, p_value, _, _ = chi2_contingency(table)
    cramers_v = float(np.sqrt(chi_square / (len(pooled) * min(table.shape[0] - 1, table.shape[1] - 1))))
    args.output_dir.mkdir(parents=True)
    _write(args.output_dir / "cohort_vs_pooled_rank3.tsv", comparison)
    receipt = {
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "comparison_rank": 3,
        "pooled_sample_count": len(pooled),
        "study_by_state_table": table.tolist(),
        "study_order": study_names,
        "state_order": state_names,
        "study_state_chi_square": float(chi_square),
        "study_state_p_value": float(p_value),
        "study_state_cramers_v": cramers_v,
        "claim_limit": "Technical concordance only; state labels are aligned post hoc and are not Barnes subtype identities.",
    }
    (args.output_dir / "comparison_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
