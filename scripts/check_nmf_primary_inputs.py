#!/usr/bin/env python3
"""Gate the four harmonised RNA matrices for primary-only and pooled NMF."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path


STUDIES = ("E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568")
EXPECTED = {"E-MTAB-7223": 9, "E-MTAB-10801": 13, "E-MTAB-11000": 11, "E-MTAB-14568": 27}


def matrix_signature(path: Path) -> tuple[list[str], list[str], str]:
    digest = hashlib.sha256()
    genes: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:2] != ["gene_id", "gene_name"]:
            raise ValueError(f"unexpected matrix header: {path}")
        for row in reader:
            genes.append(row[0])
            digest.update(row[0].encode())
            digest.update(b"\n")
    return header[2:], genes, digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--matrix-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    gene_signatures = set()
    total = 0
    run_union: set[str] = set()
    ocm_union: set[str] = set()
    patient_union: set[str] = set()
    for study in STUDIES:
        matrix = args.matrix_root / study / "gene_log1p_tpm.tsv.gz"
        qc_path = args.repo / "data/processed/rna" / study / "aggregation/sample_qc.tsv"
        samples, genes, signature = matrix_signature(matrix)
        gene_signatures.add((len(genes), signature))
        with qc_path.open(encoding="utf-8", newline="") as handle:
            qc = list(csv.DictReader(handle, delimiter="\t"))
        qc_by_run = {row["run_accession"]: row for row in qc}
        if set(samples) != set(qc_by_run) or len(qc_by_run) != len(qc):
            raise ValueError(f"matrix/QC run mismatch for {study}")
        primary = [row for row in qc if row["primary_cohort_eligible"].lower() == "true"]
        if len(primary) != EXPECTED[study]:
            raise ValueError(f"{study}: expected {EXPECTED[study]} primary runs, found {len(primary)}")
        for row in primary:
            if row["sample_class"] != "tumour" or row["histotype_group"] != "HGSOC":
                raise ValueError(f"{study}: non-HGSOC tumour passed primary gate: {row['run_accession']}")
            if row["run_accession"] in run_union or row["canonical_ocm_id"] in ocm_union:
                raise ValueError(f"duplicate pooled run/OCM: {row['run_accession']} / {row['canonical_ocm_id']}")
            run_union.add(row["run_accession"])
            ocm_union.add(row["canonical_ocm_id"])
            patient_union.add(row["patient_id"])
        total += len(primary)
        rows.append({
            "study": study,
            "matrix_samples": len(samples),
            "primary_samples": len(primary),
            "genes": len(genes),
            "gene_id_sha256": signature,
        })
    if len(gene_signatures) != 1:
        raise ValueError("gene identifiers/order differ between cohort matrices")
    if total != 60 or len(run_union) != 60 or len(ocm_union) != 60 or len(patient_union) != 52:
        raise ValueError(
            f"pooled manifest mismatch: runs={len(run_union)}, OCMs={len(ocm_union)}, patients={len(patient_union)}"
        )
    receipt = {
        "status": "PASS",
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "cohorts": rows,
        "pooled_primary_runs": 60,
        "pooled_primary_ocms": 60,
        "pooled_patients": 52,
        "sample_policy": "primary_cohort_eligible=true; tumour; HGSOC",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
