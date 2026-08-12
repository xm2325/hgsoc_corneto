"""Fail-closed preparation of the public GSE189955 RNA matrices.

The GEO cell metadata does not contain a per-cell malignant flag.  This
adapter therefore preserves the authors' reported cell classes and creates a
separate, explicitly qualified comparison role.  In particular, an epithelial
cell from an HGSOC tumour site is a *candidate* malignant epithelial cell; it
is not promoted to a definitive malignant call by this module.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
import urllib.request
import zipfile
from collections import defaultdict
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import IO, Any

from hgsoc_corneto.xlsx import read_first_worksheet

REQUIRED_METADATA_COLUMNS = {
    "",
    "orig.ident",
    "nCount_RNA",
    "nFeature_RNA",
    "location",
    "Cell_Type",
}

# These meanings are recoverable by matching the location codes in the GEO
# cell identifiers to the anatomical sampling descriptions in published Table
# S1.  They are not free-text guesses from downstream biology.
LOCATION_MAP = {
    "PT": ("primary_tumour", "primary"),
    "AS": ("ascites", "effusion"),
    "GO": ("omental_metastasis", "metastasis"),
    "LN": ("lymph_node_metastasis", "metastasis"),
    "ME": ("mesenteric_metastasis", "metastasis"),
    "PE": ("peritoneal_metastasis", "metastasis"),
    "RE": ("rectal_metastasis", "metastasis"),
    "RL": ("round_ligament_metastasis", "metastasis"),
    "FT": ("fallopian_tube", "reference"),
}

EXPECTED_CELL_TYPES = {
    "Ciliated cells",
    "Epithelial cells",
    "Fibroblasts",
    "Macrophages",
    "Secretory cells",
    "T cells",
}


def _open_text(path: Path, mode: str) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str) -> str:
    observed = sha256(path)
    if observed != expected:
        raise ValueError(
            f"SHA-256 mismatch for {path}: expected {expected}, observed {observed}"
        )
    return observed


def download(url: str, target: Path) -> None:
    """Atomically download one public source file."""

    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as handle:
        temporary = Path(handle.name)
        request = urllib.request.Request(url, headers={"User-Agent": "hgsoc-corneto/0.1"})
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                shutil.copyfileobj(response, handle)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(target)


def fetch_sources(config: dict[str, Any], raw_dir: Path) -> dict[str, Path]:
    """Download and checksum the two GEO RNA files plus published Table S1."""

    raw_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for key in ("metadata", "counts"):
        spec = config[key]
        path = raw_dir / spec["filename"]
        if not path.exists():
            download(spec["url"], path)
        verify_sha256(path, spec["sha256"])
        paths[key] = path

    supplement = config["supplement"]
    archive_path = raw_dir / supplement["filename"]
    if not archive_path.exists():
        download(supplement["url"], archive_path)
    table_path = raw_dir / supplement["member"]
    with zipfile.ZipFile(archive_path) as archive:
        members = archive.namelist()
        if members.count(supplement["member"]) != 1:
            raise ValueError(
                f"Expected exactly one {supplement['member']} in supplement; found {members}"
            )
        with archive.open(supplement["member"]) as source, tempfile.NamedTemporaryFile(
            dir=raw_dir, delete=False
        ) as handle:
            temporary = Path(handle.name)
            shutil.copyfileobj(source, handle)
        temporary.replace(table_path)
    verify_sha256(table_path, supplement["member_sha256"])
    paths["supplement_archive"] = archive_path
    paths["patient_table"] = table_path
    return paths


def parse_patient_table_rows(rows: list[list[object]]) -> list[dict[str, object]]:
    """Parse the simple, published Table S1 grid into patient-level records."""

    header_index = next(
        (
            index
            for index, row in enumerate(rows)
            if row and str(row[0]).strip() == "Patient ID"
        ),
        None,
    )
    if header_index is None:
        raise ValueError("Patient ID header not found in GSE189955 Table S1")

    output: list[dict[str, object]] = []
    for row in rows[header_index + 2 :]:
        if not row or row[0] in (None, ""):
            continue
        patient_id = str(row[0]).strip()
        if patient_id == "Sum" or not re.fullmatch(r"OC\d{2}", patient_id):
            continue
        if len(row) < 12:
            raise ValueError(f"Incomplete Table S1 row for {patient_id}: {row}")
        disease = str(row[2]).strip()
        anatomical_samples = str(row[11]).strip()
        output.append(
            {
                "patient_id": patient_id,
                "age_at_diagnosis": int(float(row[1])),
                "disease_reported": disease,
                "stage_reported": str(row[3]).strip(),
                "anatomical_samples_reported": anatomical_samples,
                "is_hgsoc": disease == "HGSC",
                "is_fallopian_tube_reference_donor": "Fallopian tube" in anatomical_samples,
                "source": "PMC9627134 Table S1",
            }
        )
    if len({row["patient_id"] for row in output}) != len(output):
        raise ValueError("Duplicate patient IDs in GSE189955 Table S1")
    return output


def parse_patient_table(path: Path) -> list[dict[str, object]]:
    return parse_patient_table_rows(read_first_worksheet(path))


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def _count_matrix_cell_id(metadata_cell_id: str) -> str:
    """Reconcile two documented naming differences between the GEO files.

    The processed metadata spells the negative EpCAM gate as ``Epcam_n`` and
    seven mouth-sort identifiers as ``mouth_n2`` ... ``mouth_n8``.  The count
    matrix spells the same identifiers ``Epcam-`` and ``mouth-2`` ...
    ``mouth-8``.  No fuzzy matching is performed beyond these exact patterns.
    """

    output = metadata_cell_id.replace("_Epcam_n_", "_Epcam-_")
    return re.sub(r"_mouth_n([2-8])$", r"_mouth-\1", output)


def _comparison_role(
    *, is_hgsoc: bool, is_ft_reference: bool, site_code: str, cell_type: str
) -> tuple[str, str, str]:
    epithelial = cell_type in {"Epithelial cells", "Secretory cells", "Ciliated cells"}
    if is_hgsoc and site_code != "FT" and cell_type == "Epithelial cells":
        return (
            "hgsoc_epithelial_candidate",
            "candidate",
            "No per-cell malignant flag is present in GEO RNA metadata; HGSOC disease, "
            "non-FT site, and the published epithelial annotation are indirect evidence.",
        )
    if is_ft_reference and site_code == "FT" and epithelial:
        return (
            "normal_ft_epithelial_reference",
            "reference",
            "The six donors had no fallopian-tube malignancy, but several had other "
            "gynaecological diseases; this is not a healthy-population reference.",
        )
    if is_hgsoc and site_code != "FT" and cell_type == "Fibroblasts":
        return (
            "hgsoc_site_fibroblast_reference",
            "reference",
            "Published metadata report fibroblasts, not a validated CAF subtype; interpret "
            "this group as a stromal/CAF-proxy comparison only.",
        )
    return (
        "context_only",
        "not_assigned",
        "Not one of the pre-specified HGSOC epithelial, FT epithelial, or fibroblast contrasts.",
    )


def load_cell_metadata(
    path: Path, patient_rows: list[dict[str, object]]
) -> list[dict[str, object]]:
    """Load source cell metadata and attach qualified patient/site annotations."""

    patients = {str(row["patient_id"]): row for row in patient_rows}
    output: list[dict[str, object]] = []
    seen: set[str] = set()
    with _open_text(path, "rt") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or []) != REQUIRED_METADATA_COLUMNS:
            raise ValueError(
                "Unexpected GSE189955 metadata columns: "
                f"expected {sorted(REQUIRED_METADATA_COLUMNS)}, got {reader.fieldnames}"
            )
        for source in reader:
            cell_id = source[""].strip()
            patient_id = source["orig.ident"].strip()
            site_code = source["location"].strip()
            cell_type = source["Cell_Type"].strip()
            if not cell_id or cell_id in seen:
                raise ValueError(f"Blank or duplicate cell ID: {cell_id!r}")
            seen.add(cell_id)
            if patient_id not in patients:
                raise ValueError(f"Cell {cell_id} has patient absent from Table S1: {patient_id}")
            if site_code not in LOCATION_MAP:
                raise ValueError(f"Unknown location code for {cell_id}: {site_code}")
            if cell_type not in EXPECTED_CELL_TYPES:
                raise ValueError(f"Unknown cell type for {cell_id}: {cell_type}")
            if int(source["nCount_RNA"]) <= 0 or int(source["nFeature_RNA"]) <= 0:
                raise ValueError(f"Non-positive source QC value for {cell_id}")
            patient = patients[patient_id]
            role, malignancy_status, claim_limit = _comparison_role(
                is_hgsoc=bool(patient["is_hgsoc"]),
                is_ft_reference=bool(patient["is_fallopian_tube_reference_donor"]),
                site_code=site_code,
                cell_type=cell_type,
            )
            site_label, site_category = LOCATION_MAP[site_code]
            group_id = "__".join(
                ["GSE189955", patient_id, site_label, _slug(cell_type)]
            )
            output.append(
                {
                    "cell_id": cell_id,
                    "count_matrix_cell_id": _count_matrix_cell_id(cell_id),
                    "cell_id_reconciliation": (
                        "exact"
                        if _count_matrix_cell_id(cell_id) == cell_id
                        else "documented_geo_file_spelling_difference"
                    ),
                    "patient_id": patient_id,
                    "n_count_rna_reported": int(source["nCount_RNA"]),
                    "n_feature_rna_reported": int(source["nFeature_RNA"]),
                    "site_code": site_code,
                    "site_label": site_label,
                    "site_category": site_category,
                    "cell_type_reported": cell_type,
                    "disease_reported": patient["disease_reported"],
                    "stage_reported": patient["stage_reported"],
                    "is_hgsoc": patient["is_hgsoc"],
                    "is_fallopian_tube_reference_donor": patient[
                        "is_fallopian_tube_reference_donor"
                    ],
                    "comparison_role": role,
                    "malignancy_status": malignancy_status,
                    "claim_limit": claim_limit,
                    "pseudobulk_group_id": group_id,
                }
            )
    return output


def build_group_metadata(cells: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for cell in cells:
        groups[str(cell["pseudobulk_group_id"])].append(cell)
    output: list[dict[str, object]] = []
    for group_id in sorted(groups):
        members = groups[group_id]
        invariant_fields = [
            "patient_id",
            "site_code",
            "site_label",
            "site_category",
            "cell_type_reported",
            "disease_reported",
            "stage_reported",
            "is_hgsoc",
            "is_fallopian_tube_reference_donor",
            "comparison_role",
            "malignancy_status",
            "claim_limit",
        ]
        record: dict[str, object] = {"pseudobulk_group_id": group_id}
        for field in invariant_fields:
            values = {member[field] for member in members}
            if len(values) != 1:
                raise ValueError(f"Non-invariant {field} in {group_id}: {values}")
            record[field] = next(iter(values))
        record["n_cells"] = len(members)
        record["sum_n_count_rna_reported"] = sum(
            int(member["n_count_rna_reported"]) for member in members
        )
        output.append(record)
    return output


def write_tsv(rows: list[dict[str, object]], path: Path) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def aggregate_pseudobulk(
    counts_path: Path,
    cells: list[dict[str, object]],
    output_path: Path,
    *,
    expected_genes: int | None = None,
) -> dict[str, int]:
    """Stream the gene-by-cell matrix and sum raw counts by patient/site/type."""

    group_ids = sorted({str(cell["pseudobulk_group_id"]) for cell in cells})
    group_index = {group_id: index for index, group_id in enumerate(group_ids)}
    matrix_to_cell: dict[str, dict[str, object]] = {}
    for cell in cells:
        matrix_id = str(cell["count_matrix_cell_id"])
        if matrix_id in matrix_to_cell:
            raise ValueError(f"Non-unique reconciled count-matrix cell ID: {matrix_id}")
        matrix_to_cell[matrix_id] = cell
    gene_count = 0
    input_count_sum = 0
    output_count_sum = 0
    seen_genes: set[str] = set()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(counts_path, "rt") as source, _open_text(output_path, "wt") as target:
        reader = csv.reader(source)
        writer = csv.writer(target)
        header = next(reader, None)
        if header is None or header[0] not in ("", "gene"):
            raise ValueError("Count matrix lacks the expected leading gene column")
        matrix_cell_ids = header[1:]
        if len(matrix_cell_ids) != len(set(matrix_cell_ids)):
            raise ValueError("Duplicate cell IDs in count-matrix header")
        if set(matrix_cell_ids) != set(matrix_to_cell):
            metadata_only = sorted(set(matrix_to_cell) - set(matrix_cell_ids))[:10]
            counts_only = sorted(set(matrix_cell_ids) - set(matrix_to_cell))[:10]
            raise ValueError(
                "Count-matrix cells do not exactly match reconciled metadata cells; "
                f"metadata-only={metadata_only}, counts-only={counts_only}"
            )
        ordered_cells = [matrix_to_cell[cell_id] for cell_id in matrix_cell_ids]
        cell_group_indices = [
            group_index[str(cell["pseudobulk_group_id"])] for cell in ordered_cells
        ]
        writer.writerow(["gene", *group_ids])
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"Count row {line_number} has {len(row)} fields; expected {len(header)}"
                )
            gene = row[0].strip()
            if not gene or gene in seen_genes:
                raise ValueError(f"Blank or duplicate gene at row {line_number}: {gene!r}")
            seen_genes.add(gene)
            sums = [0] * len(group_ids)
            for raw_count, index in zip(row[1:], cell_group_indices, strict=True):
                try:
                    count = int(raw_count)
                except ValueError as error:
                    raise ValueError(
                        f"Non-integer count for {gene} at row {line_number}: {raw_count!r}"
                    ) from error
                if count < 0:
                    raise ValueError(f"Negative count for {gene} at row {line_number}")
                sums[index] += count
                input_count_sum += count
            writer.writerow([gene, *sums])
            output_count_sum += sum(sums)
            gene_count += 1
    if expected_genes is not None and gene_count != expected_genes:
        raise ValueError(f"Expected {expected_genes} genes, observed {gene_count}")
    if input_count_sum != output_count_sum:
        raise ValueError("Pseudobulk aggregation did not preserve the total count sum")
    return {
        "cells": len(cells),
        "groups": len(group_ids),
        "genes": gene_count,
        "input_count_sum": input_count_sum,
        "output_count_sum": output_count_sum,
    }


def audit_expected_counts(
    config: dict[str, Any],
    patients: list[dict[str, object]],
    cells: list[dict[str, object]],
) -> None:
    expected = config["expected"]
    observed = {
        "cells": len(cells),
        "patients": len(patients),
        "hgsoc_patients": sum(bool(row["is_hgsoc"]) for row in patients),
        "fallopian_tube_reference_donors": sum(
            bool(row["is_fallopian_tube_reference_donor"]) for row in patients
        ),
    }
    for key, value in observed.items():
        if value != int(expected[key]):
            raise ValueError(f"Expected {expected[key]} {key}, observed {value}")


def prepare(
    config: dict[str, Any], paths: dict[str, Path], output_dir: Path
) -> dict[str, object]:
    patients = parse_patient_table(paths["patient_table"])
    cells = load_cell_metadata(paths["metadata"], patients)
    audit_expected_counts(config, patients, cells)
    groups = build_group_metadata(cells)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(patients, output_dir / "patient_metadata.tsv")
    write_tsv(cells, output_dir / "cell_metadata_audited.tsv")
    write_tsv(groups, output_dir / "pseudobulk_sample_metadata.tsv")
    aggregation = aggregate_pseudobulk(
        paths["counts"],
        cells,
        output_dir / "patient_celltype_site_raw_counts.csv.gz",
        expected_genes=int(config["expected"]["genes"]),
    )
    receipt: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "dataset": config["accession"],
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "sources": {
            key: {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for key, path in sorted(paths.items())
        },
        "aggregation": aggregation,
        "patients": len(patients),
        "hgsoc_patients": sum(bool(row["is_hgsoc"]) for row in patients),
        "fallopian_tube_reference_donors": sum(
            bool(row["is_fallopian_tube_reference_donor"]) for row in patients
        ),
        "comparison_group_counts": {
            role: sum(row["comparison_role"] == role for row in groups)
            for role in sorted({str(row["comparison_role"]) for row in groups})
        },
        "cell_id_reconciliation": {
            status: sum(cell["cell_id_reconciliation"] == status for cell in cells)
            for status in sorted({str(cell["cell_id_reconciliation"]) for cell in cells})
        },
        "claim_limits": [
            "GEO RNA metadata do not provide a per-cell malignant flag.",
            "Published fibroblast labels are a stromal/CAF proxy, not a CAF subtype call.",
            "Fallopian-tube references are from donors without FT malignancy, not an "
            "unselected healthy population.",
            "The GEO Series lists 20 GSM records while its processed RNA metadata contain "
            "22 patient identities; all 22 are validated against published Table S1.",
            "Multiple regions of the same patient and site cannot be separated because the "
            "processed RNA metadata expose only a location code, not a region identifier.",
            "Cells are technical observations; patient-level pseudobulks are the units for "
            "inference.",
        ],
    }
    receipt_path = output_dir / "receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def load_config(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def prepare_from_existing_files(
    config_path: Path, raw_dir: Path, output_dir: Path
) -> dict[str, object]:
    config = load_config(config_path)
    paths = {
        "metadata": raw_dir / config["metadata"]["filename"],
        "counts": raw_dir / config["counts"]["filename"],
        "supplement_archive": raw_dir / config["supplement"]["filename"],
        "patient_table": raw_dir / config["supplement"]["member"],
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required GSE189955 source files: {missing}")
    verify_sha256(paths["metadata"], config["metadata"]["sha256"])
    verify_sha256(paths["counts"], config["counts"]["sha256"])
    verify_sha256(paths["patient_table"], config["supplement"]["member_sha256"])
    return prepare(config, paths, output_dir)


def patient_ids(rows: Iterable[dict[str, object]]) -> set[str]:
    """Small public helper used by source audits."""

    return {str(row["patient_id"]) for row in rows}
