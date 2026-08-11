"""MAGE-TAB and ENA run-report harmonisation."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any

from .identifiers import apply_alias, parse_source_name

STUDY_PROJECTS = {
    "E-MTAB-7223": "PRJEB28709",
    "E-MTAB-10801": "PRJEB46736",
    "E-MTAB-11000": "PRJEB47842",
    "E-MTAB-14568": "PRJEB81794",
}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _first(row: dict[str, str], keys: list[str]) -> str | None:
    for key in keys:
        value = row.get(key, "").strip()
        if value:
            return value
    return None


def _sample_class(source_name: str, parsed_material: str | None, row: dict[str, str]) -> str:
    if source_name.upper() in {"FNE", "KURAMOCHI"}:
        return "cell_line_control"
    if parsed_material in {"tumour", "stroma"}:
        return parsed_material
    context = " ".join(
        value.lower()
        for key, value in row.items()
        if value
        and (
            "cell type" in key.lower()
            or "sampling site" in key.lower()
            or "source name" == key.lower()
        )
    )
    if "stromal" in context or "fibroblast" in context:
        return "stroma"
    if "tumour" in context or "tumor" in context or "epithelial" in context:
        return "tumour"
    return "unknown"


def load_rna_runs(raw_root: str | Path, aliases: dict[str, str]) -> list[dict[str, Any]]:
    """Return one validated record per ENA run across all four accessions."""

    root = Path(raw_root)
    output: list[dict[str, Any]] = []
    for study, project in STUDY_PROJECTS.items():
        sdrf_path = root / "biostudies" / f"{study}.sdrf.txt"
        ena_path = root / "ena" / f"{project}.read_run.tsv"
        sdrf_rows = _read_tsv(sdrf_path)
        ena_rows = _read_tsv(ena_path)
        ena_by_run = {row["run_accession"]: row for row in ena_rows}
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in sdrf_rows:
            grouped[row["Comment[ENA_RUN]"].strip()].append(row)

        if set(grouped) != set(ena_by_run):
            missing_sdrf = sorted(set(ena_by_run) - set(grouped))
            missing_ena = sorted(set(grouped) - set(ena_by_run))
            raise ValueError(
                f"{study} SDRF/ENA run mismatch: SDRF-missing={missing_sdrf}; "
                f"ENA-missing={missing_ena}"
            )

        for run in sorted(grouped):
            rows = grouped[run]
            ena = ena_by_run[run]
            if len(rows) != 2:
                raise ValueError(f"{run}: expected two paired-end SDRF rows, found {len(rows)}")
            source_names = {row["Source Name"].strip() for row in rows}
            if len(source_names) != 1:
                raise ValueError(f"{run}: multiple source names {source_names}")
            source_name = source_names.pop()
            parsed = parse_source_name(source_name)
            canonical = apply_alias(parsed.source_biospecimen_id, aliases)
            sample_class = _sample_class(source_name, parsed.material, rows[0])
            fastq_ftp = [value for value in ena["fastq_ftp"].split(";") if value]
            fastq_md5 = [value for value in ena["fastq_md5"].split(";") if value]
            fastq_bytes = [int(value) for value in ena["fastq_bytes"].split(";") if value]
            if not (len(fastq_ftp) == len(fastq_md5) == len(fastq_bytes) == 2):
                raise ValueError(f"{run}: incomplete paired FASTQ metadata")
            if ena["library_layout"] != "PAIRED":
                raise ValueError(f"{run}: unexpected layout {ena['library_layout']}")

            output.append(
                {
                    "study_accession": study,
                    "ena_project": project,
                    "run_accession": run,
                    "experiment_accession": ena["experiment_accession"],
                    "ena_sample_accession": ena["sample_accession"],
                    "secondary_sample_accession": ena["secondary_sample_accession"],
                    "biosd_sample_accession": rows[0].get("Comment[BioSD_SAMPLE]", ""),
                    "source_name": source_name,
                    "ena_sample_alias": ena["sample_alias"],
                    "ena_sample_title": ena["sample_title"],
                    "source_biospecimen_id": parsed.source_biospecimen_id,
                    "canonical_ocm_id": canonical,
                    "patient_id": parsed.patient_id,
                    "sample_class": sample_class,
                    "passage": parsed.passage,
                    "individual_reported": _first(
                        rows[0],
                        ["Characteristics[individual]", "Factor Value[individual]"],
                    ),
                    "disease_reported": _first(
                        rows[0],
                        ["Characteristics[disease]", "Factor Value[disease]"],
                    ),
                    "organism_part_reported": _first(
                        rows[0],
                        ["Characteristics[organism part]", "Characteristics[sampling site]"],
                    ),
                    "cell_type_reported": _first(
                        rows[0],
                        ["Characteristics[cell type]", "Factor Value[cell type]"],
                    ),
                    "library_strategy": ena["library_strategy"],
                    "library_selection": ena["library_selection"],
                    "library_layout": ena["library_layout"],
                    "instrument_model": ena["instrument_model"],
                    "read_count": int(ena["read_count"]),
                    "base_count": int(ena["base_count"]),
                    "fastq_ftp": ";".join(fastq_ftp),
                    "fastq_md5": ";".join(fastq_md5),
                    "fastq_bytes": ";".join(str(value) for value in fastq_bytes),
                    "fastq_total_bytes": sum(fastq_bytes),
                    "first_public": ena["first_public"],
                    "metadata_source": f"BioStudies {study} SDRF + ENA {project} read_run",
                }
            )

    if len(output) != 117:
        raise ValueError(f"Expected 117 RNA runs, found {len(output)}")
    if len({row["run_accession"] for row in output}) != 117:
        raise ValueError("Duplicate run accessions")
    return output
