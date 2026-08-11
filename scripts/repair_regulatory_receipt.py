#!/usr/bin/env python3
"""Replace errored sample rows with validated single-sample retry receipts."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


class RepairError(ValueError):
    """Raised when a retry receipt cannot safely repair the original."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RepairError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise RepairError(f"{path}: root must be an object")
    return value


def repair(original_path: Path, retry_paths: list[Path]) -> dict[str, Any]:
    original = _read(original_path)
    samples = original.get("samples")
    if not isinstance(samples, list) or not samples:
        raise RepairError("original samples missing")
    by_run: dict[str, dict[str, Any]] = {}
    for row in samples:
        if not isinstance(row, dict) or not isinstance(row.get("run_accession"), str):
            raise RepairError("original has an invalid sample row")
        run = row["run_accession"]
        if run in by_run:
            raise RepairError(f"original has duplicate sample {run}")
        by_run[run] = row
    repair_log = []
    for retry_path in retry_paths:
        retry = _read(retry_path)
        for field in ("study", "primary_only", "method", "source_sha256"):
            if retry.get(field) != original.get(field):
                raise RepairError(f"{retry_path}: {field} differs from original")
        retry_samples = retry.get("samples")
        if not isinstance(retry_samples, list) or not retry_samples:
            raise RepairError(f"{retry_path}: samples missing")
        for replacement in retry_samples:
            if not isinstance(replacement, dict):
                raise RepairError(f"{retry_path}: invalid retry sample")
            run = replacement.get("run_accession")
            if run not in by_run:
                raise RepairError(f"{retry_path}: retry sample {run!r} is absent from original")
            if by_run[run].get("status") != "error":
                raise RepairError(f"{retry_path}: original sample {run} is not errored")
            if replacement.get("status") == "error":
                raise RepairError(f"{retry_path}: replacement sample {run} is still errored")
            old = by_run[run]
            by_run[run] = copy.deepcopy(replacement)
            repair_log.append(
                {
                    "run_accession": run,
                    "original_status": old.get("status"),
                    "replacement_status": replacement.get("status"),
                    "retry_receipt": str(retry_path),
                    "retry_sha256": _sha256(retry_path),
                }
            )
    if not repair_log:
        raise RepairError("no sample was repaired")
    remaining_errors = sorted(
        run for run, row in by_run.items() if row.get("status") == "error"
    )
    if remaining_errors:
        raise RepairError(f"unrepaired sample errors remain: {remaining_errors}")
    result = copy.deepcopy(original)
    result["schema_version"] = "corneto_regulatory_pilot_repaired.v1"
    result["samples"] = [by_run[row["run_accession"]] for row in samples]
    statuses = [row.get("status") for row in result["samples"]]
    result["status"] = (
        "completed"
        if any(status in {"optimal", "optimal_inaccurate"} for status in statuses)
        else "blocked"
    )
    result["repair"] = {
        "response_blind": True,
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "original_receipt": str(original_path),
        "original_sha256": _sha256(original_path),
        "replacements": repair_log,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--retry", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = repair(args.original, args.retry)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    except (RepairError, OSError) as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": result["status"],
                "repaired_samples": len(result["repair"]["replacements"]),
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
