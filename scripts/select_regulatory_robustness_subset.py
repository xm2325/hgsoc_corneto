#!/usr/bin/env python3
"""Select a deterministic, stratified response-blind robustness subset.

The subset is a pilot for parameter stability, not the final cohort result.
It keeps prespecified longitudinal anchors when available, then balances
regulatory status and NMF state before filling remaining slots by OCM/run ID.
The original manifest is never modified.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any


ANCHORS = ("OCM231", "OCM341", "OCM66")


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"missing header: {path}")
        return list(reader.fieldnames), list(reader)


def _truth(value: str | None) -> bool:
    return str(value).casefold() in {"true", "1", "yes"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--regulatory", type=Path, required=True)
    parser.add_argument("--nmf", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--study", default="E-MTAB-14568")
    parser.add_argument("--n", type=int, default=12)
    args = parser.parse_args()
    if args.n < len(ANCHORS):
        raise SystemExit("--n must leave room for prespecified anchors")

    fields, rows = _read_tsv(args.manifest)
    rows = [
        row
        for row in rows
        if row.get("study_accession") == args.study
        and _truth(row.get("primary_cohort_eligible"))
        and row.get("sample_class") == "tumour"
    ]
    by_run = {row.get("run_accession", ""): row for row in rows}
    reg_status: dict[str, str] = {}
    try:
        payload: Any = json.loads(args.regulatory.read_text(encoding="utf-8"))
        for sample in payload.get("samples", []):
            if isinstance(sample, dict) and sample.get("run_accession"):
                reg_status[str(sample["run_accession"])] = str(sample.get("status", "unknown"))
    except (OSError, json.JSONDecodeError):
        pass

    _, nmf_rows = _read_tsv(args.nmf)
    nmf_state = {
        row.get("run_accession", ""): row.get("independent_state", "unknown")
        for row in nmf_rows
        if row.get("run_accession")
    }

    def key(row: dict[str, str]) -> tuple[str, str, str, str]:
        run = row.get("run_accession", "")
        return (
            reg_status.get(run, "unknown"),
            nmf_state.get(run, "unknown"),
            row.get("canonical_ocm_id", ""),
            run,
        )

    selected: list[dict[str, str]] = []
    selected_runs: set[str] = set()
    for anchor in ANCHORS:
        candidates = sorted(
            (row for row in rows if row.get("canonical_ocm_id") == anchor),
            key=key,
        )
        if candidates:
            selected.append(candidates[0])
            selected_runs.add(candidates[0].get("run_accession", ""))

    strata: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        run = row.get("run_accession", "")
        if run not in selected_runs:
            strata[(reg_status.get(run, "unknown"), nmf_state.get(run, "unknown"))].append(row)
    for values in strata.values():
        values.sort(key=key)
    # Round-robin strata so one state/status does not dominate the pilot.
    while len(selected) < args.n and strata:
        progressed = False
        for stratum in sorted(list(strata)):
            values = strata.get(stratum, [])
            if values:
                row = values.pop(0)
                selected.append(row)
                selected_runs.add(row.get("run_accession", ""))
                progressed = True
                if len(selected) >= args.n:
                    break
            if not values:
                strata.pop(stratum, None)
        if not progressed:
            break
    if len(selected) < args.n:
        remaining = sorted((row for row in rows if row.get("run_accession", "") not in selected_runs), key=key)
        selected.extend(remaining[: args.n - len(selected)])
    selected = selected[: args.n]
    selected.sort(key=lambda row: row.get("run_accession", ""))
    if len(selected) != args.n:
        raise SystemExit(f"only {len(selected)} eligible rows; requested {args.n}")

    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.output_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    selection = {
        "status": "selected",
        "study": args.study,
        "n": args.n,
        "policy": "anchors OCM231/OCM341/OCM66, then round-robin regulatory status × NMF state, deterministic OCM/run sort",
        "runs": [
            {
                "run_accession": row.get("run_accession"),
                "canonical_ocm_id": row.get("canonical_ocm_id"),
                "patient_id": row.get("patient_id"),
                "regulatory_status": reg_status.get(row.get("run_accession", ""), "unknown"),
                "nmf_state": nmf_state.get(row.get("run_accession", ""), "unknown"),
            }
            for row in selected
        ],
    }
    args.output_json.write_text(json.dumps(selection, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "selected", "n": args.n, "output_manifest": str(args.output_manifest)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
