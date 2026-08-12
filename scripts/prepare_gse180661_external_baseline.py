#!/usr/bin/env python3
"""Gate GSE180661 metadata and optionally build patient-level pseudobulks."""

from __future__ import annotations

import argparse
from pathlib import Path

from hgsoc_corneto.external.gse180661 import (
    aggregate_10x_h5_pseudobulk,
    load_source_manifest,
    prepare_metadata_gate,
)
from hgsoc_corneto.io import write_json

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "config/external/gse180661_sources.json",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=ROOT / "data/raw/external/GSE180661",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "data/processed/external/GSE180661",
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        required=True,
        help="Pre-specified minimum cells per patient x cell-type x site group.",
    )
    parser.add_argument(
        "--aggregate",
        action="store_true",
        help="After the metadata gate, validate the HDF5 and sum raw counts.",
    )
    parser.add_argument("--chunk-cells", type=int, default=4096)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = load_source_manifest(args.manifest)
    receipt = prepare_metadata_gate(
        manifest_path=args.manifest,
        cell_metadata_path=args.raw_dir / "GSE180661_GEO_cells.tsv.gz",
        sample_inventory_path=(
            args.raw_dir / "41586_2022_5496_MOESM3_ESM.xlsx"
        ),
        mutational_signatures_path=(
            args.raw_dir / "41586_2022_5496_MOESM4_ESM.xlsx"
        ),
        output_dir=args.output_dir,
        min_cells=args.min_cells,
    )
    if args.aggregate:
        matrix_spec = next(
            item for item in manifest["files"] if item["role"] == "count_matrix"
        )
        aggregation = aggregate_10x_h5_pseudobulk(
            matrix_path=args.raw_dir / matrix_spec["filename"],
            cell_map_path=args.output_dir / "cell_to_patient_celltype_site.tsv.gz",
            group_table_path=args.output_dir / "patient_celltype_site_groups.tsv",
            output_path=args.output_dir / "patient_celltype_site_raw_counts.tsv.gz",
            expected_bytes=int(matrix_spec["bytes"]),
            frozen_sha256=matrix_spec["sha256"],
            chunk_cells=args.chunk_cells,
        )
        write_json(args.output_dir / "pseudobulk_receipt.json", aggregation)
        print(
            f"GSE180661: status={aggregation['status']} "
            f"cells={aggregation['included_cells']} groups={aggregation['eligible_groups']}"
        )
    else:
        print(
            f"GSE180661: status={receipt['status']} "
            f"cells={receipt['observed']['cells']} "
            f"groups={receipt['observed']['pseudobulk_groups']}"
        )


if __name__ == "__main__":
    main()
