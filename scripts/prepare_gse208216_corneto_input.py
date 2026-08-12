#!/usr/bin/env python3
"""Build a gene-symbol CPM matrix and model-level manifest for GSE208216."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import tempfile
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from hgsoc_corneto.external_validation import (
    audit_gse_count_matrix,
    file_sha256,
    read_gene_map,
)
from hgsoc_corneto.io import deterministic_gzip_text_writer, write_json, write_tsv


def prepare(
    counts: Path,
    gene_map_path: Path,
    contract_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists():
        raise ValueError(f"refusing to overwrite {output_dir}")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    samples = contract["samples"]
    sample_ids = [row["sample_id"] for row in samples]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("GSE208216 contract has duplicate samples")
    audit = audit_gse_count_matrix(
        counts,
        expected_sha256=contract["count_matrix"]["sha256"],
        expected_samples=sample_ids,
        expected_gene_rows=int(contract["count_matrix"]["gene_rows"]),
    )
    raw_map = read_gene_map(gene_map_path)
    gene_map: dict[str, str] = {}
    for gene_id, symbol in raw_map.items():
        normalized = gene_id.split(".", maxsplit=1)[0]
        previous = gene_map.get(normalized)
        if previous is not None and previous != symbol:
            raise ValueError(f"conflicting symbol map for {normalized}")
        gene_map[normalized] = symbol

    symbol_counts: dict[str, list[int]] = defaultdict(lambda: [0] * len(sample_ids))
    library_sizes = [0] * len(sample_ids)
    mapped_gene_rows = 0
    with gzip.open(counts, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        matrix_samples = header[1:]
        indices = [matrix_samples.index(sample) + 1 for sample in sample_ids]
        for line, row in enumerate(reader, start=2):
            values = [int(row[index]) for index in indices]
            if any(value < 0 for value in values):
                raise ValueError(f"negative count at row {line}")
            for index, value in enumerate(values):
                library_sizes[index] += value
            symbol = gene_map.get(row[0].split(".", maxsplit=1)[0])
            if symbol:
                mapped_gene_rows += 1
                for index, value in enumerate(values):
                    symbol_counts[symbol][index] += value
    if any(size <= 0 for size in library_sizes) or not symbol_counts:
        raise ValueError("empty library or no mapped gene symbols")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output_dir.parent) as temporary_name:
        stage = Path(temporary_name)
        matrix_path = stage / "gene_symbol_cpm.tsv.gz"
        with deterministic_gzip_text_writer(matrix_path) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_name", *sample_ids])
            for symbol in sorted(symbol_counts):
                writer.writerow(
                    [
                        symbol,
                        *(
                            format(
                                symbol_counts[symbol][index]
                                * 1_000_000.0
                                / library_sizes[index],
                                ".12g",
                            )
                            for index in range(len(sample_ids))
                        ),
                    ]
                )
        manifest_rows = [
            {
                "study_accession": "GSE208216",
                "run_accession": row["sample_id"],
                "patient_id": row["sample_id"],
                "comparison_role": row["sample_class"],
                "geo_accession": row["geo_accession"],
                "analysis_unit": "organoid_model",
                "claim_limit": (
                    "Model-level public RNA profile; model count is not a patient-cohort size."
                ),
            }
            for row in samples
        ]
        manifest_path = stage / "corneto_manifest.tsv"
        write_tsv(manifest_path, manifest_rows)
        receipt = {
            "schema_version": "gse208216_corneto_input.v1",
            "status": "completed",
            "created_at_utc": datetime.now(UTC).isoformat(),
            "study_accession": "GSE208216",
            "dimensions": {
                "samples": len(sample_ids),
                "pdo_models": sum(
                    row["sample_class"] == "hgsoc_organoid" for row in samples
                ),
                "fallopian_tube_models": sum(
                    row["sample_class"] == "fallopian_tube_organoid" for row in samples
                ),
                "source_gene_rows": audit["gene_row_count"],
                "mapped_source_gene_rows": mapped_gene_rows,
                "unique_gene_symbols": len(symbol_counts),
            },
            "normalization": {
                "operation": "sum explicitly mapped Ensembl counts by symbol, then CPM",
                "scope": "GSE208216 only",
                "log_transform": "none; downstream CORNETO applies log1p",
                "library_sizes": dict(zip(sample_ids, library_sizes, strict=True)),
            },
            "sources": {
                "counts": {"path": str(counts), "sha256": file_sha256(counts)},
                "gene_map": {
                    "path": str(gene_map_path),
                    "sha256": file_sha256(gene_map_path),
                },
                "contract": {
                    "path": str(contract_path),
                    "sha256": file_sha256(contract_path),
                },
            },
            "outputs": {
                "matrix": {"filename": matrix_path.name, "sha256": file_sha256(matrix_path)},
                "manifest": {
                    "filename": manifest_path.name,
                    "sha256": file_sha256(manifest_path),
                },
            },
            "claim_limit": (
                "Eleven PDO and three FT organoid models provide in-vitro model evidence, "
                "not an independent 14-patient clinical cohort."
            ),
        }
        write_json(stage / "receipt.json", receipt)
        stage.replace(output_dir)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--gene-map", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(args.counts, args.gene_map, args.contract, args.output_dir),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
