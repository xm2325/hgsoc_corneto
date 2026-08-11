#!/usr/bin/env python3
"""Create a deterministic one-OCM-per-patient NMF sensitivity matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"run_accession", "patient_id", "canonical_ocm_id", "study_accession", "primary_cohort_eligible"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("manifest is missing patient-balancing columns")
    rows = [row for row in rows if row["primary_cohort_eligible"].strip().lower() == "true"]
    by_run = {row["run_accession"]: row for row in rows}
    if len(by_run) != len(rows):
        raise ValueError("duplicate primary run accession in manifest")
    return by_run


def prepare(matrix: Path, manifest: Path, output_dir: Path) -> dict[str, object]:
    if output_dir.exists():
        raise ValueError(f"refusing to overwrite {output_dir}")
    by_run = _read_manifest(manifest)
    by_patient: dict[str, list[dict[str, str]]] = {}
    for row in by_run.values():
        by_patient.setdefault(row["patient_id"], []).append(row)
    selected = [sorted(rows, key=lambda row: (row["study_accession"], row["run_accession"]))[0] for rows in by_patient.values()]
    selected.sort(key=lambda row: (row["study_accession"], row["run_accession"]))
    if len(selected) != 52 or len({row["patient_id"] for row in selected}) != 52:
        raise ValueError(f"expected 52 unique patients, selected {len(selected)}")
    selected_runs = {row["run_accession"] for row in selected}
    output_dir.mkdir(parents=True, exist_ok=False)
    out_matrix = output_dir / "gene_log1p_tpm_patient_balanced.tsv.gz"
    with gzip.open(matrix, "rt", encoding="utf-8", newline="") as source, gzip.open(out_matrix, "wt", encoding="utf-8", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        header = next(reader)
        if header[:2] != ["gene_id", "gene_name"]:
            raise ValueError("unexpected matrix header")
        indices = [header.index(row["run_accession"]) for row in selected]
        writer.writerow(["gene_id", "gene_name", *(header[index] for index in indices)])
        for row in reader:
            if len(row) != len(header):
                raise ValueError("matrix row width mismatch")
            writer.writerow([row[0], row[1], *(row[index] for index in indices)])
    qc_path = output_dir / "sample_qc.tsv"
    qc_fields = ["study_accession", "original_study_accession", "run_accession", "canonical_ocm_id", "patient_id", "sample_class", "histotype_group", "primary_cohort_eligible"]
    with qc_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=qc_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in selected:
            values = {field: row.get(field, "true" if field == "primary_cohort_eligible" else "") for field in qc_fields}
            # The NMF runner treats the cohort label as the QC study accession;
            # preserve the real source study separately for provenance.
            values["study_accession"] = "patient_balanced_primary"
            values["original_study_accession"] = row["study_accession"]
            writer.writerow(values)
    manifest_path = output_dir / "selected_samples.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["study_accession", "run_accession", "canonical_ocm_id", "patient_id"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in writer.fieldnames} for row in selected)
    receipt = {
        "status": "completed",
        "schema_version": "patient_balanced_nmf_input.v1",
        "selection_policy": "one primary OCM per patient; lexicographically first (study_accession, run_accession)",
        "source_matrix": {"path": str(matrix), "sha256": _sha(matrix)},
        "source_manifest": {"path": str(manifest), "sha256": _sha(manifest)},
        "selected_run_count": len(selected),
        "selected_patient_count": len({row["patient_id"] for row in selected}),
        "selected_study_counts": {study: sum(row["study_accession"] == study for row in selected) for study in sorted({row["study_accession"] for row in selected})},
        "output_matrix": str(out_matrix),
        "output_qc": str(qc_path),
        "claim_limit": "patient-balanced technical NMF sensitivity; no phenotype or causal interpretation",
    }
    (output_dir / "patient_balanced_input_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.matrix, args.manifest, args.output_dir), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
