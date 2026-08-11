#!/usr/bin/env python3
"""Build source-derived Tighe and RNA metadata tables."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hgsoc_corneto.identifiers import read_aliases  # noqa: E402
from hgsoc_corneto.io import checksum_rows, write_json, write_tsv  # noqa: E402
from hgsoc_corneto.magetab import load_rna_runs  # noqa: E402
from hgsoc_corneto.manifest import build_master_manifest  # noqa: E402
from hgsoc_corneto.tighe import parse_tighe_table_s1  # noqa: E402
from hgsoc_corneto.xlsx import parse_tighe_abcb1  # noqa: E402

TIGHE_FIELDS = [
    "table_row",
    "patient_numeric",
    "patient_id",
    "figo_stage",
    "histotype_reported",
    "histotype_group",
    "primary_tp53",
    "ocm_id_reported",
    "canonical_ocm_id",
    "chemo_naive_at_biopsy",
    "biopsy_type",
    "ocm_tp53_dna",
    "ocm_tp53_protein",
    "p53_without_nutlin3",
    "p53_with_nutlin3",
    "references",
    "source",
    "source_pmc",
    "source_pdf_page",
]

ABCB1_FIELDS = [
    "ocm_id_reported",
    "canonical_ocm_id",
    "abcb1_normalized_read_count",
    "abcb1_missing_reason",
    "source",
    "source_pmc",
]

RNA_FIELDS = [
    "study_accession",
    "ena_project",
    "run_accession",
    "experiment_accession",
    "ena_sample_accession",
    "secondary_sample_accession",
    "biosd_sample_accession",
    "source_name",
    "ena_sample_alias",
    "ena_sample_title",
    "source_biospecimen_id",
    "canonical_ocm_id",
    "patient_id",
    "sample_class",
    "passage",
    "individual_reported",
    "disease_reported",
    "organism_part_reported",
    "cell_type_reported",
    "library_strategy",
    "library_selection",
    "library_layout",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "fastq_total_bytes",
    "first_public",
    "metadata_source",
]

MASTER_ADDITIONAL_FIELDS = [
    "in_tighe_83_ocm_screen",
    "histotype_reported",
    "histotype_group",
    "figo_stage",
    "chemo_naive_at_biopsy",
    "biopsy_type",
    "tighe_table_row",
    "abcb1_normalized_read_count",
    "abcb1_missing_reason",
    "is_representative_rna_library",
    "hgsoc_tumour_eligible",
    "primary_cohort_eligible",
    "exact_paclitaxel_auc_available",
    "exact_paclitaxel_gi50_available",
    "exact_cumulative_paclitaxel_exposure_available",
    "paclitaxel_auc",
    "paclitaxel_gi50_nm",
    "cumulative_paclitaxel_exposure_mg",
]


def _typed_tighe_from_tsv(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for source in csv.DictReader(handle, delimiter="\t"):
            row: dict[str, Any] = dict(source)
            value = source.get("chemo_naive_at_biopsy", "").lower()
            row["chemo_naive_at_biopsy"] = (
                True if value == "true" else False if value == "false" else None
            )
            rows.append(row)
    return rows


def _longitudinal_rows(tighe_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[str]] = defaultdict(list)
    for row in tighe_rows:
        grouped[row["patient_id"]].append(row["canonical_ocm_id"])

    explicit_types = {
        "OCM64": "mixed_longitudinal_and_same_biopsy_fraction",
        "OCM66": "longitudinal",
        "OCM74": "longitudinal",
        "OCM110": "longitudinal",
        "OCM118": "longitudinal",
        "OCM124": "longitudinal",
        "OCM231": "longitudinal",
        "OCM288": "longitudinal",
        "OCM296": "longitudinal",
        "OCM327": "longitudinal",
        "OCM333": "longitudinal",
        "OCM341": "longitudinal",
        "OCM361": "spatial_same_time_not_longitudinal",
    }
    anchors = {
        "OCM66": "stable_control",
        "OCM231": "acquired_resistance_anchor",
        "OCM341": "acquired_resistance_anchor",
    }
    output = []
    for patient_id, models in sorted(grouped.items(), key=lambda item: int(item[0][3:])):
        if len(models) < 2:
            continue
        output.append(
            {
                "patient_id": patient_id,
                "ocm_ids": ";".join(models),
                "n_ocms": len(models),
                "relationship_type": explicit_types.get(
                    patient_id, "repeated_patient_timing_not_yet_curated"
                ),
                "analysis_role": anchors.get(patient_id, "none_prespecified"),
                "source": "Tighe2025 Table S1; relationship labels from frozen contract",
            }
        )
    return output


def _audit(
    rna_runs: list[dict[str, Any]],
    tighe_rows: list[dict[str, Any]],
    abcb1_rows: list[dict[str, Any]],
    manifest: list[dict[str, Any]],
) -> dict[str, Any]:
    by_study: dict[str, dict[str, int]] = {}
    for study in sorted({row["study_accession"] for row in rna_runs}):
        rows = [row for row in rna_runs if row["study_accession"] == study]
        by_study[study] = {
            "runs": len(rows),
            "fastq_files": sum(len(row["fastq_ftp"].split(";")) for row in rows),
            "fastq_bytes": sum(row["fastq_total_bytes"] for row in rows),
        }
    unresolved = sorted(
        {
            row["canonical_ocm_id"]
            for row in manifest
            if row["sample_class"] == "tumour"
            and row["canonical_ocm_id"]
            and not row["in_tighe_83_ocm_screen"]
        }
    )
    duplicate_groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in manifest:
        if row["canonical_ocm_id"]:
            duplicate_groups[(row["canonical_ocm_id"], row["sample_class"])].append(
                row["source_name"]
            )
    duplicates = {
        f"{ocm_id}|{sample_class}": sorted(source_names)
        for (ocm_id, sample_class), source_names in duplicate_groups.items()
        if len(source_names) > 1
    }
    return {
        "rna": {
            "studies": by_study,
            "runs_total": len(rna_runs),
            "fastq_files_total": sum(len(row["fastq_ftp"].split(";")) for row in rna_runs),
            "fastq_bytes_total": sum(row["fastq_total_bytes"] for row in rna_runs),
            "sample_class_counts": dict(Counter(row["sample_class"] for row in rna_runs)),
        },
        "tighe": {
            "table_s1_ocms": len(tighe_rows),
            "table_s1_patients": len({row["patient_id"] for row in tighe_rows}),
            "table_s2_ocms": len(abcb1_rows),
            "hgsoc_ocms": sum(row["histotype_group"] == "HGSOC" for row in tighe_rows),
        },
        "manifest": {
            "rows": len(manifest),
            "primary_cohort_rna_runs": sum(row["primary_cohort_eligible"] for row in manifest),
            "primary_cohort_ocms": len(
                {row["canonical_ocm_id"] for row in manifest if row["primary_cohort_eligible"]}
            ),
            "tumour_ocms_without_tighe_table_s1_match": unresolved,
            "canonical_ocms_with_multiple_rna_runs": duplicates,
        },
        "phenotype_boundary": {
            "exact_paclitaxel_auc_rows": 0,
            "exact_paclitaxel_gi50_rows": 0,
            "exact_cumulative_paclitaxel_exposure_rows": 0,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, default=ROOT / "data/raw/metadata")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/processed/metadata")
    parser.add_argument("--aliases", type=Path, default=ROOT / "config/ocm_aliases.tsv")
    parser.add_argument("--tighe-pdf", type=Path, default=ROOT / "tmp/pdfs/mmc1.pdf")
    parser.add_argument(
        "--tighe-table",
        type=Path,
        help="Use a previously extracted Table S1 TSV instead of parsing the PDF",
    )
    parser.add_argument(
        "--tighe-abcb1",
        type=Path,
        default=ROOT / "data/raw/metadata/mmc2.xlsx",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    aliases = read_aliases(args.aliases)
    tighe_rows = (
        _typed_tighe_from_tsv(args.tighe_table)
        if args.tighe_table
        else parse_tighe_table_s1(args.tighe_pdf)
    )
    abcb1_rows = parse_tighe_abcb1(args.tighe_abcb1)
    rna_runs = load_rna_runs(args.raw_root, aliases)
    manifest = build_master_manifest(rna_runs, tighe_rows, abcb1_rows)

    write_tsv(args.output_dir / "tighe_ocm_characteristics.tsv", tighe_rows, TIGHE_FIELDS)
    write_tsv(args.output_dir / "tighe_abcb1.tsv", abcb1_rows, ABCB1_FIELDS)
    write_tsv(args.output_dir / "rna_runs.tsv", rna_runs, RNA_FIELDS)
    write_tsv(
        args.output_dir / "ocm_master_manifest.tsv",
        manifest,
        RNA_FIELDS + MASTER_ADDITIONAL_FIELDS,
    )
    write_tsv(
        args.output_dir / "longitudinal_families.tsv",
        _longitudinal_rows(tighe_rows),
    )
    write_json(
        args.output_dir / "manifest_audit.json",
        _audit(rna_runs, tighe_rows, abcb1_rows, manifest),
    )
    raw_files = [path for path in args.raw_root.rglob("*") if path.is_file()]
    write_tsv(
        args.output_dir / "raw_metadata_checksums.tsv",
        checksum_rows(raw_files, ROOT),
        ["path", "bytes", "sha256"],
    )


if __name__ == "__main__":
    main()
