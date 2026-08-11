#!/usr/bin/env python3
"""Build a provenance-checked pooled expression matrix for the primary OCM cohort.

The four study matrices contain controls, stroma, non-HGSOC samples, and alternate
libraries.  This gate validates every matrix against the frozen master manifest,
then streams only ``primary_cohort_eligible`` columns into one deterministic gzip
matrix.  No normalization or batch correction is performed here.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import tempfile
from collections import Counter
from contextlib import ExitStack
from datetime import UTC, datetime
from pathlib import Path
from typing import TextIO

from hgsoc_corneto.io import (
    deterministic_gzip_text_writer,
    read_tsv,
    sha256,
    write_json,
    write_tsv,
)


def _parse_matrix_spec(value: str) -> tuple[str, Path]:
    study, separator, raw_path = value.partition("=")
    if not separator or not study or not raw_path:
        raise argparse.ArgumentTypeError("matrix must be STUDY=PATH")
    return study, Path(raw_path)


def _open_matrix(stack: ExitStack, path: Path) -> TextIO:
    if path.suffix == ".gz":
        return stack.enter_context(gzip.open(path, "rt", encoding="utf-8", newline=""))
    return stack.enter_context(path.open(encoding="utf-8", newline=""))


def _is_true(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def _validate_primary(rows: list[dict[str, str]]) -> None:
    required = {
        "study_accession",
        "run_accession",
        "canonical_ocm_id",
        "patient_id",
        "sample_class",
        "histotype_group",
        "is_representative_rna_library",
        "primary_cohort_eligible",
    }
    missing = required - set(rows[0]) if rows else required
    if missing:
        raise ValueError(f"manifest missing required fields: {sorted(missing)}")
    for row in rows:
        if row["sample_class"] != "tumour":
            raise ValueError(f"primary run is not tumour: {row['run_accession']}")
        if row["histotype_group"] != "HGSOC":
            raise ValueError(f"primary run is not HGSOC: {row['run_accession']}")
        if not _is_true(row["is_representative_rna_library"]):
            raise ValueError(f"primary run is not representative: {row['run_accession']}")
        for field in ("run_accession", "canonical_ocm_id", "patient_id"):
            if not row[field] or row[field] == "NA":
                raise ValueError(f"primary run has missing {field}: {row['run_accession']}")


def build_pooled_matrix(
    *,
    manifest_path: Path,
    matrices: list[tuple[str, Path]],
    output_dir: Path,
    value_name: str,
    expected_samples: int,
    expected_patients: int,
) -> dict[str, object]:
    if output_dir.exists():
        raise ValueError(f"refusing to overwrite existing output: {output_dir}")
    if len({study for study, _ in matrices}) != len(matrices):
        raise ValueError("duplicate study matrix arguments")
    manifest = read_tsv(manifest_path)
    manifest_by_study: dict[str, list[dict[str, str]]] = {}
    seen_runs: set[str] = set()
    for row in manifest:
        run = row["run_accession"]
        if run in seen_runs:
            raise ValueError(f"duplicate run in manifest: {run}")
        seen_runs.add(run)
        manifest_by_study.setdefault(row["study_accession"], []).append(row)
    requested_studies = [study for study, _ in matrices]
    unknown = set(requested_studies) - set(manifest_by_study)
    if unknown:
        raise ValueError(f"matrix studies absent from manifest: {sorted(unknown)}")
    primary = [
        row
        for study in requested_studies
        for row in manifest_by_study[study]
        if _is_true(row["primary_cohort_eligible"])
    ]
    _validate_primary(primary)
    if len(primary) != expected_samples:
        raise ValueError(f"expected {expected_samples} primary runs, found {len(primary)}")
    patient_count = len({row["patient_id"] for row in primary})
    if patient_count != expected_patients:
        raise ValueError(f"expected {expected_patients} patients, found {patient_count}")
    ocm_count = len({row["canonical_ocm_id"] for row in primary})
    if ocm_count != expected_samples:
        raise ValueError(
            f"primary run-to-OCM mapping is not one-to-one: {len(primary)} runs, {ocm_count} OCMs"
        )
    primary_by_run = {row["run_accession"]: row for row in primary}

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    output_matrix = stage / f"gene_{value_name}.tsv.gz"
    output_manifest = stage / "sample_manifest.tsv"
    try:
        with ExitStack() as stack:
            readers: list[tuple[str, Path, csv.reader, list[str], list[int]]] = []
            ordered_primary: list[dict[str, str]] = []
            for study, path in matrices:
                if not path.is_file():
                    raise FileNotFoundError(path)
                reader = csv.reader(_open_matrix(stack, path), delimiter="\t")
                try:
                    header = next(reader)
                except StopIteration as error:
                    raise ValueError(f"empty matrix: {path}") from error
                if header[:2] != ["gene_id", "gene_name"]:
                    raise ValueError(f"unexpected matrix header: {path}")
                samples = header[2:]
                if len(samples) != len(set(samples)):
                    raise ValueError(f"duplicate sample columns: {path}")
                expected_runs = {row["run_accession"] for row in manifest_by_study[study]}
                if set(samples) != expected_runs:
                    missing = sorted(expected_runs - set(samples))
                    extra = sorted(set(samples) - expected_runs)
                    raise ValueError(
                        f"{study} matrix/manifest mismatch; missing={missing[:5]}, extra={extra[:5]}"
                    )
                indices = [index for index, run in enumerate(samples) if run in primary_by_run]
                selected_rows = [primary_by_run[samples[index]] for index in indices]
                ordered_primary.extend(selected_rows)
                readers.append((study, path, reader, header, indices))

            if len({row["run_accession"] for row in ordered_primary}) != expected_samples:
                raise ValueError("pooled primary run set is incomplete or duplicated")
            with deterministic_gzip_text_writer(output_matrix) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["gene_id", "gene_name", *(row["run_accession"] for row in ordered_primary)]
                )
                gene_count = 0
                while True:
                    rows: list[list[str] | None] = []
                    for _, _, reader, _, _ in readers:
                        rows.append(next(reader, None))
                    if all(row is None for row in rows):
                        break
                    if any(row is None for row in rows):
                        raise ValueError("study matrices have different gene-row counts")
                    assert all(row is not None for row in rows)
                    typed_rows = [row for row in rows if row is not None]
                    gene_key = typed_rows[0][:2]
                    if any(row[:2] != gene_key for row in typed_rows[1:]):
                        raise ValueError(f"study matrices have different gene order at row {gene_count + 2}")
                    values: list[str] = []
                    for row, (_, path, _, header, indices) in zip(typed_rows, readers, strict=True):
                        if len(row) != len(header):
                            raise ValueError(f"wrong field count in {path} at row {gene_count + 2}")
                        values.extend(row[index + 2] for index in indices)
                    writer.writerow([*gene_key, *values])
                    gene_count += 1
        manifest_fields = [
            "study_accession",
            "run_accession",
            "canonical_ocm_id",
            "patient_id",
            "sample_class",
            "histotype_group",
            "chemo_naive_at_biopsy",
            "biopsy_type",
            "primary_cohort_eligible",
        ]
        write_tsv(output_manifest, ordered_primary, manifest_fields)
        by_study = Counter(row["study_accession"] for row in ordered_primary)
        receipt: dict[str, object] = {
            "status": "completed",
            "completed_at_utc": datetime.now(UTC).isoformat(),
            "value_name": value_name,
            "normalization": "none; primary columns copied verbatim from study matrices",
            "gene_alignment": "exact gene_id and gene_name order required across all matrices",
            "gene_count": gene_count,
            "primary_run_count": len(ordered_primary),
            "primary_ocm_count": ocm_count,
            "patient_count": patient_count,
            "study_counts": dict(sorted(by_study.items())),
            "eligibility": "primary_cohort_eligible=true; tumour; HGSOC; representative RNA library",
            "inputs": {
                "manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
                "matrices": [
                    {"study_accession": study, "path": str(path), "sha256": sha256(path)}
                    for study, path in matrices
                ],
            },
            "outputs": {
                "matrix": {
                    "path": str(output_dir / output_matrix.name),
                    "bytes": output_matrix.stat().st_size,
                    "sha256": sha256(output_matrix),
                },
                "sample_manifest": {
                    "path": str(output_dir / output_manifest.name),
                    "sha256": sha256(output_manifest),
                },
            },
            "claim_limit": "pooled technical input only; no batch correction, subtype, or response claim",
        }
        write_json(stage / "pooled_expression_receipt.json", receipt)
        os.replace(stage, output_dir)
        return receipt
    except Exception:
        raise RuntimeError(f"pooled-expression staging retained for inspection: {stage}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--matrix", action="append", type=_parse_matrix_spec, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--value-name", choices=("tpm", "log1p_tpm"), required=True)
    parser.add_argument("--expected-samples", type=int, default=60)
    parser.add_argument("--expected-patients", type=int, default=52)
    args = parser.parse_args()
    receipt = build_pooled_matrix(
        manifest_path=args.manifest,
        matrices=args.matrix,
        output_dir=args.output_dir,
        value_name=args.value_name,
        expected_samples=args.expected_samples,
        expected_patients=args.expected_patients,
    )
    print(json.dumps({key: receipt[key] for key in ("status", "gene_count", "primary_run_count", "patient_count", "study_counts")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
