#!/usr/bin/env python3
"""Read-only gate for Salmon quantification completeness.

Compares manifest runs with per-run quantification receipts and required Salmon
outputs. It never opens FASTQ files and never modifies quantification outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


REQUIRED = ("quant.sf", "cmd_info.json", "aux_info/meta_info.json")


def _row(study: str, run: str, root: Path) -> dict[str, object]:
    out = root / study / run
    receipt_path = out / "run_receipt.json"
    result: dict[str, object] = {
        "study_accession": study,
        "run_accession": run,
        "status": "missing",
        "quant_dir": str(out),
        "missing": [],
        "invalid": [],
    }
    if not out.is_dir():
        result["missing"] = ["quant_dir"]
        return result
    missing = [name for name in REQUIRED if not (out / name).is_file()]
    if not receipt_path.is_file():
        missing.append("run_receipt.json")
    if missing:
        result["missing"] = missing
    if receipt_path.is_file():
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            if not isinstance(receipt, dict):
                raise ValueError("receipt_not_object")
            if receipt.get("status") != "completed":
                result["invalid"].append("receipt_status")
            run_meta = receipt.get("run")
            if not isinstance(run_meta, dict):
                result["invalid"].append("receipt_run_missing")
            else:
                if run_meta.get("study_accession") != study:
                    result["invalid"].append("receipt_study_mismatch")
                if run_meta.get("run_accession") != run:
                    result["invalid"].append("receipt_run_mismatch")
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            result["invalid"].append(f"receipt_json:{type(exc).__name__}")
    if result["invalid"]:
        result["status"] = "invalid"
    elif result["missing"]:
        result["status"] = "missing"
    else:
        result["status"] = "complete"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--quant-root", type=Path, required=True)
    parser.add_argument("--study", action="append", help="Restrict to one or more studies")
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--tsv-output", type=Path)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            study, run = record["study_accession"], record["run_accession"]
            if args.study and study not in args.study:
                continue
            key = (study, run)
            if key in seen:
                raise ValueError(f"duplicate manifest run: {study}/{run}")
            seen.add(key)
            rows.append(_row(study, run, args.quant_root))
    rows.sort(key=lambda item: (str(item["study_accession"]), str(item["run_accession"])))
    counts = {status: sum(item["status"] == status for item in rows) for status in ("complete", "missing", "invalid")}
    payload = {"quant_root": str(args.quant_root), "total": len(rows), "counts": counts, "runs": rows}
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.tsv_output:
        args.tsv_output.parent.mkdir(parents=True, exist_ok=True)
        with args.tsv_output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["study_accession", "run_accession", "status", "missing", "invalid", "quant_dir"])
            for item in rows:
                writer.writerow([item["study_accession"], item["run_accession"], item["status"], ";".join(item["missing"]), ";".join(item["invalid"]), item["quant_dir"]])
    print(json.dumps({"total": len(rows), "counts": counts}, sort_keys=True))


if __name__ == "__main__":
    main()
