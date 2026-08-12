#!/usr/bin/env python3
"""Fetch and checksum the public GSE180661 work-package inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from hgsoc_corneto.external.gse180661 import fetch_sources

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "config/external/gse180661_sources.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "data/raw/external/GSE180661",
    )
    parser.add_argument(
        "--include-matrix",
        action="store_true",
        help="Also fetch the 30.3-GiB HDF5 matrix (off by default).",
    )
    parser.add_argument(
        "--write-frozen-manifest",
        type=Path,
        help=(
            "After an opt-in matrix download, write a new manifest with the "
            "observed matrix SHA-256 frozen. The source manifest is never overwritten."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    receipt = fetch_sources(
        args.manifest,
        args.output_dir,
        include_matrix=args.include_matrix,
    )
    if args.write_frozen_manifest is not None:
        if not args.include_matrix:
            raise ValueError("--write-frozen-manifest requires --include-matrix")
        from hgsoc_corneto.external.gse180661 import write_frozen_matrix_manifest

        write_frozen_matrix_manifest(
            args.manifest,
            args.write_frozen_manifest,
            receipt,
        )
    print(
        f"{receipt['study_accession']}: status={receipt['status']} "
        f"files={len(receipt['files'])} include_matrix={receipt['include_matrix']}"
    )


if __name__ == "__main__":
    main()
