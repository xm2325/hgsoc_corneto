#!/usr/bin/env python3
"""Fetch, verify, and patient-pseudobulk the public GSE189955 RNA data."""

from __future__ import annotations

import argparse
from pathlib import Path

from hgsoc_corneto.external.gse189955 import (
    fetch_sources,
    load_config,
    prepare,
    prepare_from_existing_files,
)

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "config/external/gse189955_sources.json",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=ROOT / "data/raw/external/GSE189955",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "data/processed/external/GSE189955",
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="Download missing official source files before preparation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.fetch:
        config = load_config(args.config)
        paths = fetch_sources(config, args.raw_dir)
        receipt = prepare(config, paths, args.output_dir)
    else:
        receipt = prepare_from_existing_files(args.config, args.raw_dir, args.output_dir)
    aggregation = receipt["aggregation"]
    print(
        f"{receipt['dataset']}: status={receipt['status']} "
        f"cells={aggregation['cells']} groups={aggregation['groups']} "
        f"genes={aggregation['genes']}"
    )


if __name__ == "__main__":
    main()
