#!/usr/bin/env python3
"""Build paired metadata and CORNETO-ready GSE277107 expression matrices."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hgsoc_corneto.external.gse277107 import prepare_dataset

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=ROOT / "config/external/gse277107_sources.json",
    )
    parser.add_argument(
        "--source-dir", type=Path, default=ROOT / "data/external/GSE277107/raw"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "data/external/GSE277107/processed"
    )
    args = parser.parse_args()
    receipt = prepare_dataset(
        source_manifest_path=args.source_manifest,
        source_dir=args.source_dir,
        output_dir=args.output_dir,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
