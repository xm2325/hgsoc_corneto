#!/usr/bin/env python3
"""Convert a gene-by-pseudobulk raw-count TSV to dataset-internal CPM.

The operation is intentionally separate from CORNETO so its input checksum,
library sizes and output checksum are independently auditable.  The CORNETO
loader applies log1p after reading this non-negative CPM matrix.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from hgsoc_corneto.io import deterministic_gzip_text_writer, sha256


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def scan_counts(path: Path) -> tuple[list[str], list[float], int]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 3 or header[0] not in {"gene", "gene_name"}:
            raise ValueError("raw-count matrix must be TSV with gene plus >=2 samples")
        samples = header[1:]
        if len(samples) != len(set(samples)):
            raise ValueError("raw-count matrix has duplicate samples")
        totals = [0.0] * len(samples)
        genes: set[str] = set()
        rows = 0
        for line, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"row {line} has {len(row)} fields, expected {len(header)}")
            gene = row[0].strip()
            if not gene or gene in genes:
                raise ValueError(f"blank or duplicate gene at row {line}")
            genes.add(gene)
            for index, raw in enumerate(row[1:]):
                try:
                    value = float(raw)
                except ValueError as error:
                    raise ValueError(f"non-numeric count at row {line}") from error
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f"invalid count at row {line}")
                totals[index] += value
            rows += 1
    if rows == 0 or any(value <= 0 for value in totals):
        raise ValueError("matrix is empty or contains a zero-size library")
    return samples, totals, rows


def normalize(
    input_path: Path,
    output_path: Path,
    receipt_path: Path,
    dataset: str,
    expected_samples: int | None,
    expected_genes: int | None,
) -> dict[str, object]:
    if output_path.exists() or receipt_path.exists():
        raise ValueError("refusing to overwrite CPM output or receipt")
    samples, totals, genes = scan_counts(input_path)
    if expected_samples is not None and len(samples) != expected_samples:
        raise ValueError(f"expected {expected_samples} samples, observed {len(samples)}")
    if expected_genes is not None and genes != expected_genes:
        raise ValueError(f"expected {expected_genes} genes, observed {genes}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=output_path.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with open_text(input_path) as source, deterministic_gzip_text_writer(temporary) as target:
            reader = csv.reader(source, delimiter="\t")
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            next(reader)
            writer.writerow(["gene_name", *samples])
            for row in reader:
                writer.writerow(
                    [
                        row[0],
                        *(
                            format(float(value) * 1_000_000.0 / totals[index], ".12g")
                            for index, value in enumerate(row[1:])
                        ),
                    ]
                )
        temporary.replace(output_path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    receipt = {
        "schema_version": "external_pseudobulk_cpm.v1",
        "status": "completed",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "dataset": dataset,
        "input": {"path": str(input_path), "sha256": sha256(input_path)},
        "output": {"path": str(output_path), "sha256": sha256(output_path)},
        "dimensions": {"genes": genes, "samples": len(samples)},
        "normalization": {
            "operation": "count / sample_library_sum * 1e6",
            "scope": "within this external dataset only",
            "log_transform": "none; downstream CORNETO applies log1p",
            "library_sums": dict(zip(samples, totals, strict=True)),
        },
        "claim_limit": "CPM corrects library size but not cell composition or batch effects.",
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=receipt_path.parent, delete=False
    ) as handle:
        temporary_receipt = Path(handle.name)
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    temporary_receipt.replace(receipt_path)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--expected-samples", type=int)
    parser.add_argument("--expected-genes", type=int)
    args = parser.parse_args()
    print(
        json.dumps(
            normalize(
                args.input,
                args.output,
                args.receipt,
                args.dataset,
                args.expected_samples,
                args.expected_genes,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
