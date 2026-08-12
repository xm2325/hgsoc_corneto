#!/usr/bin/env python3
"""Fetch and checksum the frozen public GSE277107 inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from hgsoc_corneto.external.gse277107 import fetch_sources

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=ROOT / "config/external/gse277107_sources.json",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "data/external/GSE277107/raw"
    )
    args = parser.parse_args()
    for role, path in sorted(fetch_sources(args.source_manifest, args.output_dir).items()):
        print(f"{role}\t{path}")


if __name__ == "__main__":
    main()
