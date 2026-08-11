#!/usr/bin/env python3
"""Aggregate a complete Salmon study from transcripts to GENCODE genes."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import tempfile
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

from hgsoc_corneto.io import (
    deterministic_gzip_text_writer,
    read_tsv,
    sha256,
    write_json,
    write_tsv,
)
from hgsoc_corneto.rna import (
    GeneRecord,
    SalmonGeneSample,
    aggregate_salmon_quant,
    iter_gene_matrix_rows,
    load_gencode_gene_map,
    load_rna_run_specs,
)


def _git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return result.stdout.strip()


def _write_matrix(
    path: Path,
    genes: tuple[GeneRecord, ...],
    samples: tuple[SalmonGeneSample, ...],
    value: str,
) -> None:
    with deterministic_gzip_text_writer(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_name", *(sample.run_accession for sample in samples)])
        for gene_id, gene_name, values in iter_gene_matrix_rows(genes, samples, value=value):
            writer.writerow([gene_id, gene_name, *(format(item, ".12g") for item in values)])


def _load_receipt(path: Path, study: str, run_accession: str) -> dict[str, object]:
    receipt = json.loads(path.read_text(encoding="utf-8"))
    run = receipt.get("run")
    if receipt.get("status") != "completed" or not isinstance(run, dict):
        raise ValueError(f"Incomplete Salmon receipt: {path}")
    if run.get("study_accession") != study or run.get("run_accession") != run_accession:
        raise ValueError(f"Salmon receipt/run mismatch: {path}")
    return receipt


def _sample_metadata(manifest: Path, study: str) -> dict[str, dict[str, str]]:
    rows = [row for row in read_tsv(manifest) if row["study_accession"] == study]
    indexed = {row["run_accession"]: row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"Duplicate run rows in {manifest}")
    return indexed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-manifest", type=Path, required=True)
    parser.add_argument("--master-manifest", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--quant-root", type=Path, required=True)
    parser.add_argument("--matrix-output-dir", type=Path, required=True)
    parser.add_argument("--report-output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.matrix_output_dir.exists() or args.report_output_dir.exists():
        raise ValueError("Existing aggregation output requires inspection; refusing to overwrite")
    runs = load_rna_run_specs(args.run_manifest, study_accession=args.study)
    master = _sample_metadata(args.master_manifest, args.study)
    if set(master) != {run.run_accession for run in runs}:
        raise ValueError("Run and master manifests disagree for selected study")
    genes, transcript_map = load_gencode_gene_map(args.gtf)

    samples: list[SalmonGeneSample] = []
    receipts: dict[str, dict[str, object]] = {}
    quant_paths: dict[str, Path] = {}
    for run in runs:
        run_dir = args.quant_root / args.study / run.run_accession
        quant_path = run_dir / "quant.sf"
        receipt_path = run_dir / "run_receipt.json"
        if not quant_path.is_file() or not receipt_path.is_file():
            raise FileNotFoundError(f"Missing completed Salmon output: {run.run_accession}")
        receipts[run.run_accession] = _load_receipt(
            receipt_path, args.study, run.run_accession
        )
        quant_paths[run.run_accession] = quant_path
        samples.append(
            aggregate_salmon_quant(
                quant_path,
                run_accession=run.run_accession,
                transcript_to_gene_index=transcript_map,
                gene_count=len(genes),
            )
        )
    sample_tuple = tuple(samples)
    unmapped_union = sorted(
        {item for sample in sample_tuple for item in sample.unmapped_transcript_ids}
    )
    if unmapped_union:
        raise ValueError(
            f"{len(unmapped_union)} Salmon transcript IDs do not map to the frozen GTF; "
            f"first={unmapped_union[:5]}"
        )

    args.matrix_output_dir.parent.mkdir(parents=True, exist_ok=True)
    args.report_output_dir.parent.mkdir(parents=True, exist_ok=True)
    matrix_stage = Path(
        tempfile.mkdtemp(prefix=f".{args.study}.aggregate.", dir=args.matrix_output_dir.parent)
    )
    report_stage = Path(
        tempfile.mkdtemp(prefix=f".{args.study}.report.", dir=args.report_output_dir.parent)
    )
    try:
        matrix_files = {
            "gene_counts": matrix_stage / "gene_counts.tsv.gz",
            "gene_tpm": matrix_stage / "gene_tpm.tsv.gz",
            "gene_log1p_tpm": matrix_stage / "gene_log1p_tpm.tsv.gz",
        }
        for value, path in (
            ("counts", matrix_files["gene_counts"]),
            ("tpm", matrix_files["gene_tpm"]),
            ("log1p_tpm", matrix_files["gene_log1p_tpm"]),
        ):
            _write_matrix(path, genes, sample_tuple, value)

        write_tsv(
            matrix_stage / "gene_metadata.tsv",
            [asdict(gene) for gene in genes],
        )
        qc_rows = []
        for sample in sample_tuple:
            metadata = master[sample.run_accession]
            receipt = receipts[sample.run_accession]
            quant_path = quant_paths[sample.run_accession]
            qc_rows.append(
                {
                    "study_accession": args.study,
                    "run_accession": sample.run_accession,
                    "canonical_ocm_id": metadata["canonical_ocm_id"],
                    "patient_id": metadata["patient_id"],
                    "sample_class": metadata["sample_class"],
                    "histotype_group": metadata["histotype_group"],
                    "primary_cohort_eligible": metadata["primary_cohort_eligible"],
                    "library_types_detected": ";".join(
                        str(value)
                        for value in (receipt.get("library_types_detected") or [])
                    ),
                    "num_processed": receipt.get("num_processed"),
                    "num_mapped": receipt.get("num_mapped"),
                    "percent_mapped": receipt.get("percent_mapped"),
                    "transcript_rows": sample.transcript_rows,
                    "mapped_transcript_rows": sample.mapped_transcript_rows,
                    "estimated_count_sum": format(sample.estimated_count_sum, ".12g"),
                    "tpm_sum": format(sample.tpm_sum, ".12g"),
                    "quant_sf_sha256": sha256(quant_path),
                }
            )
        write_tsv(report_stage / "sample_qc.tsv", qc_rows)

        matrix_receipts = {
            name: {
                "path": str(args.matrix_output_dir / path.name),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for name, path in matrix_files.items()
        }
        gene_metadata_path = matrix_stage / "gene_metadata.tsv"
        matrix_receipts["gene_metadata"] = {
            "path": str(args.matrix_output_dir / gene_metadata_path.name),
            "bytes": gene_metadata_path.stat().st_size,
            "sha256": sha256(gene_metadata_path),
        }
        receipt = {
            "status": "completed",
            "completed_at_utc": datetime.now(UTC).isoformat(),
            "repo_commit": _git_commit(),
            "study_accession": args.study,
            "run_count": len(runs),
            "gene_count": len(genes),
            "transcript_count": len(transcript_map),
            "gtf": {
                "path": str(args.gtf),
                "bytes": args.gtf.stat().st_size,
                "sha256": sha256(args.gtf),
            },
            "matrices": matrix_receipts,
            "sample_qc": {
                "path": str(args.report_output_dir / "sample_qc.tsv"),
                "sha256": sha256(report_stage / "sample_qc.tsv"),
            },
            "aggregation": "sum Salmon NumReads and TPM across versioned GENCODE transcripts",
        }
        write_json(report_stage / "aggregation_receipt.json", receipt)
        os.replace(matrix_stage, args.matrix_output_dir)
        os.replace(report_stage, args.report_output_dir)
    except Exception as error:
        raise RuntimeError(
            f"Aggregation staging retained: matrix={matrix_stage}; report={report_stage}"
        ) from error
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
