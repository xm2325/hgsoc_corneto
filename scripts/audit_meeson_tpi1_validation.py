#!/usr/bin/env python3
"""Describe the public Meeson TPI1 validation table without re-analysing it.

This audit is deliberately separate from the HGSOC reconstructions.  The table
contains reaction-level experimental/control values and a published simulation,
not matched OCM perturbations.  It reports transparent descriptive agreement
only, so it cannot be mistaken for a validation of a newly fitted OCM model.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def _number(value: str, *, row: int, field: str) -> float:
    try:
        parsed = float(value)
    except ValueError as error:
        raise ValueError(f"row {row}: {field} is not numeric: {value!r}") from error
    if not math.isfinite(parsed):
        raise ValueError(f"row {row}: {field} is non-finite")
    return parsed


def _pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2:
        return None
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    left_ss = sum((x - left_mean) ** 2 for x in left)
    right_ss = sum((y - right_mean) ** 2 for y in right)
    if left_ss == 0 or right_ss == 0:
        return None
    return numerator / math.sqrt(left_ss * right_ss)


def _summary(values: list[float]) -> dict[str, float | int]:
    return {
        "count": len(values),
        "nonzero_count": sum(value != 0 for value in values),
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
    }


def audit(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = ["reactions", "NTsiRNA", "siTPI1", "siTPI1_simulation"]
        if reader.fieldnames != required:
            raise ValueError(f"unexpected columns: {reader.fieldnames!r}")
        rows = list(reader)
    if not rows:
        raise ValueError("TPI1 validation table has no rows")

    reactions = [str(row["reactions"]) for row in rows]
    if len(reactions) != len(set(reactions)) or any(not value for value in reactions):
        raise ValueError("reaction identifiers are empty or non-unique")
    control = [_number(row["NTsiRNA"], row=index + 2, field="NTsiRNA") for index, row in enumerate(rows)]
    observed = [_number(row["siTPI1"], row=index + 2, field="siTPI1") for index, row in enumerate(rows)]
    simulated = [
        _number(row["siTPI1_simulation"], row=index + 2, field="siTPI1_simulation")
        for index, row in enumerate(rows)
    ]
    changed_observed = [index for index, (a, b) in enumerate(zip(control, observed)) if a != b]
    changed_simulated = [index for index, (a, b) in enumerate(zip(control, simulated)) if a != b]
    agreement = [
        index
        for index in range(len(rows))
        if (observed[index] - control[index]) == (simulated[index] - control[index])
    ]
    return {
        "status": "completed",
        "response_blind": True,
        "input": {"path": str(path), "columns": required, "reaction_count": len(rows)},
        "channels": {
            "NTsiRNA": _summary(control),
            "siTPI1": _summary(observed),
            "siTPI1_simulation": _summary(simulated),
        },
        "comparison": {
            "observed_vs_simulation_pearson": _pearson(observed, simulated),
            "control_vs_observed_pearson": _pearson(control, observed),
            "control_vs_simulation_pearson": _pearson(control, simulated),
            "observed_changed_reaction_count": len(changed_observed),
            "simulated_changed_reaction_count": len(changed_simulated),
            "exact_delta_agreement_count": len(agreement),
            "exact_delta_agreement_fraction": len(agreement) / len(rows),
            "changed_observed_reaction_ids": [reactions[index] for index in changed_observed],
            "changed_simulated_reaction_ids": [reactions[index] for index in changed_simulated],
        },
        "claim_limit": (
            "This is an audit of Meeson's published reaction-level table. It is not an "
            "independent TPI1 perturbation experiment and does not validate any new HGSOC "
            "OCM-specific reconstruction."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"status": result["status"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
