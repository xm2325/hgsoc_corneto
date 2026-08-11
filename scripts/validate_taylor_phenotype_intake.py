#!/usr/bin/env python3
"""Run the fail-closed Taylor phenotype intake gate; never run association."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

from hgsoc_corneto.phenotype import (
    PhenotypeIntakeError,
    blocked_receipt,
    validate_phenotype_intake,
)


def _write_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as out:
        temporary = Path(out.name)
        json.dump(value, out, indent=2, sort_keys=True)
        out.write("\n")
    os.replace(temporary, path)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phenotype",
        type=Path,
        default=root / "data/raw/taylor/taylor_exact_phenotypes.tsv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "data/processed/metadata/ocm_master_manifest.tsv",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=root / "config/taylor_phenotype_intake.json",
    )
    parser.add_argument(
        "--receipt",
        type=Path,
        default=root / "data/processed/taylor/phenotype_intake_gate.json",
    )
    args = parser.parse_args()

    if not args.phenotype.is_file():
        receipt = blocked_receipt(
            phenotype_path=args.phenotype,
            manifest_path=args.manifest,
            schema_path=args.schema,
        )
    else:
        try:
            receipt = validate_phenotype_intake(
                phenotype_path=args.phenotype,
                manifest_path=args.manifest,
                schema_path=args.schema,
            )
        except PhenotypeIntakeError as error:
            receipt = {
                "schema_version": 1,
                "status": "invalid_phenotype_intake",
                "association_allowed": False,
                "reason": str(error),
                "inputs": {
                    "phenotype": str(args.phenotype),
                    "manifest": str(args.manifest),
                    "schema": str(args.schema),
                },
                "phenotype_values_written": False,
                "association_run": False,
            }
    _write_json(args.receipt, receipt)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["status"] != "ready":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
