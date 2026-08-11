#!/usr/bin/env python3
"""Filter the frozen 60-condition regulatory bundle to one OCM per patient."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare(bundle: Path, selected: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise ValueError(f"refusing to overwrite {output}")
    root = json.loads(bundle.read_text(encoding="utf-8"))
    if root.get("status") != "completed" or root.get("schema_version") != "regulatory_multisample_input.v1":
        raise ValueError("unsupported input bundle")
    with selected.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    runs = [row.get("run_accession", "") for row in rows]
    if len(runs) != 52 or len(set(runs)) != 52:
        raise ValueError(f"expected 52 unique selected runs, found {len(runs)}")
    by_run = {str(row["run_accession"]): row for row in root.get("conditions", [])}
    if set(runs) - set(by_run):
        raise ValueError("selected run missing from regulatory bundle")
    conditions = [by_run[run] for run in runs]
    result = dict(root)
    result["analysis_mode"] = "patient_balanced_pooled"
    result["conditions"] = conditions
    result["input_counts"] = dict(result.get("input_counts", {}))
    result["input_counts"]["primary_conditions"] = len(conditions)
    result["input_counts"]["included_conditions"] = sum(row.get("preprocessing_status") == "included" for row in conditions)
    result["input_counts"]["blocked_conditions"] = sum(row.get("preprocessing_status") != "included" for row in conditions)
    result["input_counts"]["patient_balanced_patients"] = len({row.get("patient_id") for row in conditions})
    result["patient_balancing"] = {
        "selected_samples_path": str(selected),
        "selected_samples_sha256": _sha(selected),
        "policy": "one primary OCM per patient; inherited deterministic selection from patient-balanced NMF input",
    }
    result["source_bundle"] = {"path": str(bundle), "sha256": _sha(bundle)}
    result["claim_limit"] = "response-blind patient-balanced regulatory sensitivity; no phenotype or causal interpretation"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.bundle, args.selected, args.output)
    print(json.dumps({"status": result["status"], "included_conditions": result["input_counts"]["included_conditions"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
