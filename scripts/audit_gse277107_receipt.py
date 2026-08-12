#!/usr/bin/env python3
"""Audit a GSE277107 preparation receipt and preserve a gate result."""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

from hgsoc_corneto.external.gse277107 import audit_prepared_dataset
from hgsoc_corneto.io import write_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit_prepared_dataset(args.processed_dir)
    except Exception as exc:
        result = {
            "schema_version": "gse277107_receipt_gate.v1",
            "status": "failed",
            "audited_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
            "study_accession": "GSE277107",
            "scientific_success": False,
            "error_class": type(exc).__name__,
            "error_message": str(exc)[:1000],
        }
        write_json(args.output, result)
        raise
    write_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
