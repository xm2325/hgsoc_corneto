"""Provenance-first preparation of the public GSE180661 HGSOC atlas.

The analysis unit produced here is a patient by author-reported cell type by
anatomical site pseudobulk.  Flow-sorting gates and cells are never promoted to
independent biological replicates.  Mutational-process labels are copied from
published Supplementary Table 3, with ``FBI`` preserved as the source label and
rendered as ``foldback_inversion`` only in a separate normalized field.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import tempfile
import urllib.request
import zipfile
from collections import Counter
from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TextIO

from hgsoc_corneto.io import (
    deterministic_gzip_text_writer,
    sha256,
    write_json,
    write_tsv,
)
from hgsoc_corneto.xlsx import read_first_worksheet


class GSE180661Error(ValueError):
    """Raised when a GSE180661 source or biological-unit gate fails closed."""


CELL_METADATA_COLUMNS = (
    "cell_id",
    "sample",
    "cell_type",
    "percent.mt",
    "nCount_RNA",
    "nFeature_RNA",
    "umap50_1",
    "umap50_2",
    "cluster_label",
    "cluster_label_sub",
    "cell_type_super",
    "patient_id",
    "tumor_subsite",
    "tumor_site",
    "tumor_supersite",
    "sort_parameters",
    "therapy",
    "surgery",
)

INVENTORY_REQUIRED = {
    "spectrum_aliquot_id",
    "spectrum_sample_id",
    "patient_id",
    "cancer_type",
    "sample_type",
    "tumor_supersite",
    "tumor_site",
    "tumor_subsite",
    "tumor_type",
    "therapy",
    "procedure",
    "platform",
    "sort_parameters",
}

SIGNATURE_COLUMNS = (
    "patient_id",
    "consensus_signature",
    "consensus_hr_status",
    "wgs_signature",
    "myriad_signature",
    "myriad_gis_score",
)

SIGNATURE_MAP = {
    "HRD-Dup": "hrd_dup",
    "HRD-Del": "hrd_del",
    "FBI": "foldback_inversion",
    "Undetermined": "undetermined",
}

CELL_TYPE_TO_SUPER = {
    "Ovarian.cancer.cell": "Ovarian.cancer.super",
    "T.cell": "T.super",
    "Myeloid.cell": "Myeloid.super",
    "Fibroblast": "Fibroblast.super",
    "Plasma.cell": "B.super",
    "B.cell": "B.super",
    "Endothelial.cell": "Endothelial.super",
    "Dendritic.cell": "Myeloid.super",
    "Mast.cell": "Myeloid.super",
    "Other": "Other.super",
}

# The released per-cell table contains author-normalized anatomical labels,
# while Supplementary Table 2 retains finer surgical specimen descriptions.
# These are the complete observed non-identity pairs in the frozen sources;
# no fuzzy or ontology-based recoding is allowed.
ALLOWED_SITE_RECODINGS = {
    "tumor_subsite": {
        ("Right Adnexa", "Right Ovary"),
        ("Pelvic Peritoneum", "Left Pelvic Peritoneum"),
        ("Left Adnexa", "Left Ovary"),
        ("Bowel", "Transverse Colon"),
        ("Right Adnexa", "Right Fallopian Tube"),
        ("Lymph Node", "Left Pararenal Lymph Node with Colonic Mesentary"),
        ("Pelvic Peritoneum", "Left Pelvic Sidewall"),
    },
    "tumor_site": {("Other", "Peritoneum")},
    "tumor_supersite": {
        ("UQ", "Upper Quadrant"),
        ("Other", "Peritoneum"),
    },
}

HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")


def load_source_manifest(path: Path) -> dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "hgsoc_external_sources.v1":
        raise GSE180661Error("unsupported source-manifest schema")
    if manifest.get("study_accession") != "GSE180661":
        raise GSE180661Error("source manifest is not for GSE180661")
    files = manifest.get("files")
    if not isinstance(files, list) or {item.get("role") for item in files} != {
        "cell_metadata",
        "count_matrix",
    }:
        raise GSE180661Error("manifest must define cell_metadata and count_matrix")
    for item in files:
        if Path(str(item.get("filename", ""))).name != item.get("filename"):
            raise GSE180661Error(f"unsafe source filename for {item.get('role')}")
        if not str(item.get("url", "")).startswith(
            "https://ftp.ncbi.nlm.nih.gov/geo/"
        ):
            raise GSE180661Error(f"non-NCBI source URL for {item.get('role')}")
        digest = item.get("sha256")
        if digest is not None and re.fullmatch(r"[0-9a-f]{64}", str(digest)) is None:
            raise GSE180661Error(f"invalid SHA-256 for {item.get('role')}")
    supplement = manifest.get("supplement", {})
    if not str(supplement.get("url", "")).startswith(
        "https://www.ebi.ac.uk/europepmc/"
    ):
        raise GSE180661Error("supplement must use the Europe PMC endpoint")
    roles = {item.get("role") for item in supplement.get("members", [])}
    if roles != {"scrna_sample_inventory", "patient_mutational_signatures"}:
        raise GSE180661Error("manifest must define both required supplement members")
    return manifest


def _download(url: str, target: Path) -> None:
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


def _verify_regular_file(path: Path, spec: dict[str, Any]) -> dict[str, Any]:
    if not path.is_file():
        raise GSE180661Error(f"missing {spec['role']} source: {path}")
    observed_bytes = path.stat().st_size
    observed_sha256 = sha256(path)
    if observed_bytes != int(spec["bytes"]):
        raise GSE180661Error(
            f"byte-size mismatch for {spec['role']}: {observed_bytes} != {spec['bytes']}"
        )
    expected_sha256 = spec.get("sha256")
    if expected_sha256 is not None and observed_sha256 != expected_sha256:
        raise GSE180661Error(
            f"SHA-256 mismatch for {spec['role']}: {observed_sha256}"
        )
    if spec["role"] == "count_matrix":
        with path.open("rb") as handle:
            if handle.read(len(HDF5_SIGNATURE)) != HDF5_SIGNATURE:
                raise GSE180661Error("count matrix does not have an HDF5 signature")
    return {
        "role": spec["role"],
        "path": str(path),
        "bytes": observed_bytes,
        "sha256": observed_sha256,
        "upstream_sha256_frozen": expected_sha256 is not None,
    }


def fetch_sources(
    manifest_path: Path,
    output_dir: Path,
    *,
    include_matrix: bool = False,
) -> dict[str, Any]:
    """Fetch public inputs atomically and emit observed-digest provenance.

    The 30.3-GiB matrix is intentionally opt-in.  GEO publishes no digest for
    it, so its observed SHA-256 must be copied into a frozen manifest revision
    before a scientific run can be considered reproducible.
    """

    manifest = load_source_manifest(manifest_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    verified: list[dict[str, Any]] = []
    for spec in manifest["files"]:
        if spec["role"] == "count_matrix" and not include_matrix:
            continue
        target = output_dir / spec["filename"]
        if not target.exists():
            _download(spec["url"], target)
        verified.append(_verify_regular_file(target, spec))

    supplement = manifest["supplement"]
    archive_path = output_dir / supplement["archive_filename"]
    if not archive_path.exists():
        _download(supplement["url"], archive_path)
    with zipfile.ZipFile(archive_path) as archive:
        names = archive.namelist()
        for member in supplement["members"]:
            filename = member["filename"]
            if names.count(filename) != 1:
                raise GSE180661Error(
                    f"expected exactly one {filename} in supplementary archive"
                )
            target = output_dir / filename
            with archive.open(filename) as source, tempfile.NamedTemporaryFile(
                dir=output_dir, delete=False
            ) as handle:
                temporary = Path(handle.name)
                shutil.copyfileobj(source, handle)
            if sha256(temporary) != member["sha256"]:
                temporary.unlink(missing_ok=True)
                raise GSE180661Error(f"SHA-256 mismatch for supplement member {filename}")
            temporary.replace(target)
            verified.append(
                {
                    "role": member["role"],
                    "path": str(target),
                    "bytes": target.stat().st_size,
                    "sha256": sha256(target),
                    "upstream_sha256_frozen": True,
                }
            )

    receipt = {
        "schema_version": "gse180661_fetch_receipt.v1",
        "study_accession": "GSE180661",
        "status": "verified",
        "created_at": datetime.now(UTC).isoformat(),
        "manifest_sha256": sha256(manifest_path),
        "include_matrix": include_matrix,
        "files": sorted(verified, key=lambda item: item["role"]),
        "matrix_identity_frozen": any(
            item["role"] == "count_matrix" and item["upstream_sha256_frozen"]
            for item in verified
        ),
    }
    write_json(output_dir / "fetch_receipt.json", receipt)
    return receipt


def write_frozen_matrix_manifest(
    source_manifest_path: Path,
    output_path: Path,
    fetch_receipt: dict[str, Any],
) -> dict[str, Any]:
    """Write a new manifest revision containing the observed H5 SHA-256.

    The reviewable source manifest is never edited in place.  This keeps the
    transition from GEO's byte-size-only identity to a locally observed,
    immutable file identity explicit in version control.
    """

    if source_manifest_path.resolve() == output_path.resolve():
        raise GSE180661Error("refusing to overwrite the source manifest in place")
    if fetch_receipt.get("status") != "verified" or not fetch_receipt.get(
        "include_matrix"
    ):
        raise GSE180661Error("a verified matrix fetch receipt is required")
    matrix_receipts = [
        item for item in fetch_receipt.get("files", []) if item.get("role") == "count_matrix"
    ]
    if len(matrix_receipts) != 1:
        raise GSE180661Error("fetch receipt must contain exactly one count matrix")
    manifest = load_source_manifest(source_manifest_path)
    matrix_spec = next(item for item in manifest["files"] if item["role"] == "count_matrix")
    matrix_spec["sha256"] = matrix_receipts[0]["sha256"]
    matrix_spec["checksum_policy"] = (
        "Locally observed SHA-256 frozen after a byte-size-gated HTTPS transfer; "
        "scientific execution requires this exact digest."
    )
    manifest["matrix_checksum_provenance"] = {
        "source": "observed after official GEO HTTPS download",
        "fetch_receipt_manifest_sha256": fetch_receipt["manifest_sha256"],
        "observed_at": fetch_receipt["created_at"],
    }
    write_json(output_path, manifest)
    return manifest


def _row_dicts(rows: list[list[object]], required: set[str]) -> list[dict[str, object]]:
    if not rows:
        raise GSE180661Error("published workbook sheet is empty")
    header = [str(value).strip() if value is not None else "" for value in rows[0]]
    if not required.issubset(header):
        raise GSE180661Error(
            f"published workbook lacks required columns: {sorted(required - set(header))}"
        )
    output: list[dict[str, object]] = []
    for values in rows[1:]:
        if not values or all(value in (None, "") for value in values):
            continue
        padded = [*values, *([None] * (len(header) - len(values)))]
        output.append(dict(zip(header, padded, strict=False)))
    return output


def parse_sample_inventory_rows(rows: list[list[object]]) -> list[dict[str, str]]:
    """Parse the scRNA sheet of published Supplementary Table 2."""

    parsed = _row_dicts(rows, INVENTORY_REQUIRED)
    output: list[dict[str, str]] = []
    for source in parsed:
        if str(source["platform"]).strip() != "10x 3' GE":
            continue
        row = {key: "" if value is None else str(value).strip() for key, value in source.items()}
        patient_id = row["patient_id"]
        aliquot_id = row["spectrum_aliquot_id"]
        if re.fullmatch(r"SPECTRUM-OV-[0-9]{3}", patient_id) is None:
            raise GSE180661Error(f"invalid patient ID in sample inventory: {patient_id}")
        if not aliquot_id.startswith(f"{patient_id}_"):
            raise GSE180661Error(f"aliquot/patient mismatch in sample inventory: {aliquot_id}")
        if row["cancer_type"] != "hgsoc" or row["therapy"] != "pre-Rx":
            raise GSE180661Error(f"unexpected disease/therapy annotation for {aliquot_id}")
        output.append(row)
    if len({row["spectrum_aliquot_id"] for row in output}) != len(output):
        raise GSE180661Error("duplicate scRNA aliquot IDs in sample inventory")
    return output


def parse_sample_inventory(path: Path) -> list[dict[str, str]]:
    return parse_sample_inventory_rows(read_first_worksheet(path))


def parse_mutational_signature_rows(rows: list[list[object]]) -> list[dict[str, str]]:
    """Parse patient-level labels from published Supplementary Table 3."""

    parsed = _row_dicts(rows, set(SIGNATURE_COLUMNS))
    output: list[dict[str, str]] = []
    for source in parsed:
        row = {
            key: "" if source.get(key) is None else str(source.get(key)).strip()
            for key in SIGNATURE_COLUMNS
        }
        signature = row["consensus_signature"]
        if signature not in SIGNATURE_MAP:
            raise GSE180661Error(f"unknown consensus mutational signature: {signature}")
        row["normalized_signature"] = SIGNATURE_MAP[signature]
        output.append(row)
    if len({row["patient_id"] for row in output}) != len(output):
        raise GSE180661Error("duplicate patient IDs in mutational-signature table")
    return output


def parse_mutational_signatures(path: Path) -> list[dict[str, str]]:
    return parse_mutational_signature_rows(read_first_worksheet(path))


def iter_annotated_cells(
    path: Path,
    inventory_rows: list[dict[str, str]],
    signature_rows: list[dict[str, str]],
) -> Iterator[dict[str, str]]:
    """Stream cells after exact sample, patient, site, sorting, and subtype joins."""

    inventory = {row["spectrum_aliquot_id"]: row for row in inventory_rows}
    signatures = {row["patient_id"]: row for row in signature_rows}
    seen: set[str] = set()
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != CELL_METADATA_COLUMNS:
            raise GSE180661Error(
                "unexpected per-cell metadata columns; source schema changed"
            )
        for source in reader:
            cell_id = source["cell_id"].strip()
            if not cell_id or cell_id in seen:
                raise GSE180661Error(f"blank or duplicate cell ID: {cell_id!r}")
            seen.add(cell_id)
            sample_id = source["sample"].strip()
            if sample_id not in inventory:
                raise GSE180661Error(f"cell sample absent from Supplementary Table 2: {sample_id}")
            sample = inventory[sample_id]
            patient_id = source["patient_id"].strip()
            if patient_id != sample["patient_id"] or patient_id not in signatures:
                raise GSE180661Error(f"patient crosswalk failure for {cell_id}")
            for source_field in ("tumor_subsite", "tumor_site", "tumor_supersite"):
                pair = (source[source_field].strip(), sample[source_field])
                if pair[0] != pair[1] and pair not in ALLOWED_SITE_RECODINGS[source_field]:
                    raise GSE180661Error(
                        f"{source_field} contradicts Supplementary Table 2 for {cell_id}"
                    )
            for source_field in ("sort_parameters", "therapy"):
                if source[source_field].strip() != sample[source_field]:
                    raise GSE180661Error(
                        f"{source_field} contradicts Supplementary Table 2 for {cell_id}"
                    )
            cell_type = source["cell_type"].strip()
            cell_type_super = source["cell_type_super"].strip()
            if CELL_TYPE_TO_SUPER.get(cell_type) != cell_type_super:
                raise GSE180661Error(
                    f"unexpected cell_type/cell_type_super pairing for {cell_id}"
                )
            signature = signatures[patient_id]
            tumor_site = source["tumor_site"].strip()
            group_id = "__".join(
                [
                    "GSE180661",
                    patient_id,
                    _slug(cell_type),
                    _slug(tumor_site),
                ]
            )
            yield {
                "cell_id": cell_id,
                "sample_id": sample_id,
                "site_specimen_id": sample["spectrum_sample_id"],
                "patient_id": patient_id,
                "cell_type_reported": cell_type,
                "cell_type_super_reported": cell_type_super,
                "tumor_subsite_cell_metadata": source["tumor_subsite"].strip(),
                "tumor_site_cell_metadata": tumor_site,
                "tumor_supersite_cell_metadata": source["tumor_supersite"].strip(),
                "tumor_subsite_inventory": sample["tumor_subsite"],
                "tumor_site_inventory": sample["tumor_site"],
                "tumor_supersite_inventory": sample["tumor_supersite"],
                "sort_parameters_reported": source["sort_parameters"].strip(),
                "therapy_reported": source["therapy"].strip(),
                "consensus_signature_reported": signature["consensus_signature"],
                "consensus_hr_status_reported": signature["consensus_hr_status"],
                "normalized_signature": signature["normalized_signature"],
                "pseudobulk_group_id": group_id,
            }


def _expected_counter(expected: dict[str, Any], key: str) -> Counter[str]:
    return Counter({str(name): int(value) for name, value in expected[key].items()})


def prepare_metadata_gate(
    *,
    manifest_path: Path,
    cell_metadata_path: Path,
    sample_inventory_path: Path,
    mutational_signatures_path: Path,
    output_dir: Path,
    min_cells: int,
) -> dict[str, Any]:
    """Build a patient × cell type × site map and fail-closed input receipt."""

    if min_cells < 1:
        raise GSE180661Error("min_cells must be a pre-specified positive integer")
    manifest = load_source_manifest(manifest_path)
    expected = manifest["expected"]
    metadata_spec = next(
        item for item in manifest["files"] if item["role"] == "cell_metadata"
    )
    _verify_regular_file(cell_metadata_path, metadata_spec)
    for role, path in (
        ("scrna_sample_inventory", sample_inventory_path),
        ("patient_mutational_signatures", mutational_signatures_path),
    ):
        member = next(item for item in manifest["supplement"]["members"] if item["role"] == role)
        if not path.is_file() or sha256(path) != member["sha256"]:
            raise GSE180661Error(f"missing or checksum-invalid {role}")

    inventory = parse_sample_inventory(sample_inventory_path)
    signatures = parse_mutational_signatures(mutational_signatures_path)
    if len(inventory) != int(expected["scrna_inventory_aliquots"]):
        raise GSE180661Error("unexpected number of scRNA aliquots in Supplementary Table 2")
    if len({row["spectrum_sample_id"] for row in inventory}) != int(
        expected["scrna_inventory_site_specimens"]
    ):
        raise GSE180661Error("unexpected number of site specimens in Supplementary Table 2")
    if len(signatures) != int(expected["mutational_signature_patients_all_modalities"]):
        raise GSE180661Error("unexpected number of patient mutational-signature rows")

    output_dir.mkdir(parents=True, exist_ok=True)
    map_path = output_dir / "cell_to_patient_celltype_site.tsv.gz"
    fields = [
        "cell_id",
        "sample_id",
        "site_specimen_id",
        "patient_id",
        "cell_type_reported",
        "cell_type_super_reported",
        "tumor_subsite_cell_metadata",
        "tumor_site_cell_metadata",
        "tumor_supersite_cell_metadata",
        "tumor_subsite_inventory",
        "tumor_site_inventory",
        "tumor_supersite_inventory",
        "sort_parameters_reported",
        "therapy_reported",
        "consensus_signature_reported",
        "consensus_hr_status_reported",
        "normalized_signature",
        "pseudobulk_group_id",
    ]
    group_cells: Counter[str] = Counter()
    cell_types: Counter[str] = Counter()
    supersites: Counter[str] = Counter()
    patient_signatures: dict[str, str] = {}
    observed_samples: set[str] = set()
    observed_specimens: set[str] = set()
    group_values: dict[str, dict[str, str]] = {}
    total_cells = 0
    with deterministic_gzip_text_writer(map_path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in iter_annotated_cells(
            cell_metadata_path, inventory, signatures
        ):
            writer.writerow(row)
            total_cells += 1
            group_id = row["pseudobulk_group_id"]
            group_cells[group_id] += 1
            cell_types[row["cell_type_reported"]] += 1
            supersites[row["tumor_supersite_cell_metadata"]] += 1
            observed_samples.add(row["sample_id"])
            observed_specimens.add(row["site_specimen_id"])
            existing_signature = patient_signatures.setdefault(
                row["patient_id"], row["consensus_signature_reported"]
            )
            if existing_signature != row["consensus_signature_reported"]:
                raise GSE180661Error("mutational signature varies within a patient")
            group_values.setdefault(group_id, row)

    if total_cells != int(expected["processed_cells"]):
        raise GSE180661Error(f"processed-cell count mismatch: {total_cells}")
    if len(observed_samples) != int(expected["processed_aliquots"]):
        raise GSE180661Error("processed-aliquot count mismatch")
    if len(observed_specimens) != int(expected["processed_site_specimens"]):
        raise GSE180661Error("processed site-specimen count mismatch")
    if len(patient_signatures) != int(expected["scrna_patients"]):
        raise GSE180661Error("processed patient count mismatch")
    if cell_types != _expected_counter(expected, "cell_type_counts"):
        raise GSE180661Error(f"cell-type census mismatch: {dict(cell_types)}")
    if supersites != _expected_counter(expected, "tumor_supersite_counts"):
        raise GSE180661Error(f"anatomical supersite census mismatch: {dict(supersites)}")
    signature_counts = Counter(patient_signatures.values())
    expected_signatures = Counter(
        {
            str(name): int(value)
            for name, value in expected["mutational_signatures_in_scrna"].items()
        }
    )
    if signature_counts != expected_signatures:
        raise GSE180661Error(f"scRNA mutational-signature census mismatch: {signature_counts}")

    groups: list[dict[str, Any]] = []
    for group_id in sorted(group_cells):
        row = group_values[group_id]
        groups.append(
            {
                "pseudobulk_group_id": group_id,
                "patient_id": row["patient_id"],
                "cell_type_reported": row["cell_type_reported"],
                "cell_type_super_reported": row["cell_type_super_reported"],
                "tumor_site_cell_metadata": row["tumor_site_cell_metadata"],
                "tumor_supersite_cell_metadata": row[
                    "tumor_supersite_cell_metadata"
                ],
                "tumor_site_inventory": row["tumor_site_inventory"],
                "tumor_supersite_inventory": row["tumor_supersite_inventory"],
                "consensus_signature_reported": row[
                    "consensus_signature_reported"
                ],
                "consensus_hr_status_reported": row[
                    "consensus_hr_status_reported"
                ],
                "normalized_signature": row["normalized_signature"],
                "n_cells": group_cells[group_id],
                "passes_min_cells": group_cells[group_id] >= min_cells,
                "statistical_unit": "patient",
            }
        )
    group_path = output_dir / "patient_celltype_site_groups.tsv"
    write_tsv(group_path, groups)

    receipt = {
        "schema_version": "gse180661_input_gate.v1",
        "study_accession": "GSE180661",
        "status": "ready_for_matrix_gate",
        "created_at": datetime.now(UTC).isoformat(),
        "source_checksums": {
            "manifest": sha256(manifest_path),
            "cell_metadata": sha256(cell_metadata_path),
            "scrna_sample_inventory": sha256(sample_inventory_path),
            "patient_mutational_signatures": sha256(mutational_signatures_path),
        },
        "observed": {
            "cells": total_cells,
            "aliquots_with_retained_cells": len(observed_samples),
            "site_specimens": len(observed_specimens),
            "patients": len(patient_signatures),
            "pseudobulk_groups": len(groups),
            "eligible_pseudobulk_groups": sum(
                bool(row["passes_min_cells"]) for row in groups
            ),
            "cell_type_counts": dict(sorted(cell_types.items())),
            "tumor_supersite_counts": dict(sorted(supersites.items())),
            "patient_mutational_signature_counts": dict(
                sorted(signature_counts.items())
            ),
        },
        "aggregation_contract": {
            "operation": "sum untransformed integer UMI counts",
            "grouping": [
                "patient_id",
                "cell_type_reported",
                "tumor_site_cell_metadata",
            ],
            "flow_sort_handling": "pool CD45+, CD45-, and unsorted cells within biological group",
            "min_cells": min_cells,
            "statistical_unit": "patient",
            "forbidden_independence_assumptions": [
                "cells as biological replicates",
                "flow-sorted aliquots as biological replicates",
                "multiple sites from one patient as independent patients",
            ],
            "primary_strata": ["HRD-Dup", "HRD-Del", "FBI"],
            "excluded_from_subtype_contrast": ["Undetermined"],
            "site_rule": (
                "compare subtypes within site or adjust/stratify by site; "
                "do not pool site effects away"
            ),
        },
        "outputs": {
            "cell_map": str(map_path),
            "group_table": str(group_path),
        },
    }
    write_json(output_dir / "input_gate_receipt.json", receipt)
    return receipt


def inspect_10x_h5_matrix(
    matrix_path: Path,
    cell_map_path: Path,
    *,
    expected_bytes: int,
    frozen_sha256: str | None,
) -> dict[str, Any]:
    """Fail-closed 10x-HDF5 schema and barcode gate.

    This function intentionally imports ``h5py`` only at runtime.  A cluster
    environment must provide it; absence is a gate failure, not permission to
    guess the matrix layout.
    """

    if matrix_path.stat().st_size != expected_bytes:
        raise GSE180661Error("count-matrix byte size does not match the source manifest")
    observed_sha256 = sha256(matrix_path)
    if frozen_sha256 is None:
        raise GSE180661Error(
            "count-matrix SHA-256 is not frozen in the source manifest; record the "
            f"observed digest before scientific execution: {observed_sha256}"
        )
    if observed_sha256 != frozen_sha256:
        raise GSE180661Error("count-matrix SHA-256 does not match the frozen manifest")
    try:
        import h5py  # type: ignore[import-not-found]
    except ImportError as exc:
        raise GSE180661Error("h5py is required for the 10x-HDF5 matrix gate") from exc

    cell_ids: set[str] = set()
    with _open_text(cell_map_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell_ids.add(row["cell_id"])
    with h5py.File(matrix_path, "r") as handle:
        if "matrix" not in handle:
            raise GSE180661Error("HDF5 source is not in the expected 10x matrix layout")
        matrix = handle["matrix"]
        required = {"data", "indices", "indptr", "shape", "barcodes", "features"}
        if not required.issubset(matrix.keys()):
            raise GSE180661Error("10x HDF5 matrix is missing required datasets")
        shape = tuple(int(value) for value in matrix["shape"][:])
        barcodes = {
            value.decode("utf-8") if isinstance(value, bytes) else str(value)
            for value in matrix["barcodes"][:]
        }
        if len(shape) != 2 or shape[1] != len(barcodes):
            raise GSE180661Error("10x matrix shape/barcode dimensions disagree")
        if barcodes != cell_ids:
            raise GSE180661Error(
                "10x HDF5 barcodes do not exactly match the checksum-gated cell map"
            )
        sample = matrix["data"][: min(100000, matrix["data"].shape[0])]
        if any(float(value) < 0 or not float(value).is_integer() for value in sample):
            raise GSE180661Error("matrix data are not untransformed non-negative counts")
    return {
        "status": "ready_for_pseudobulk_sum",
        "matrix_sha256": observed_sha256,
        "genes": shape[0],
        "cells": shape[1],
        "barcode_identity": "exact",
        "matrix_layout": "10x_hdf5_csc",
    }


def _decode_h5_strings(values: Any) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def aggregate_10x_h5_pseudobulk(
    *,
    matrix_path: Path,
    cell_map_path: Path,
    group_table_path: Path,
    output_path: Path,
    expected_bytes: int,
    frozen_sha256: str | None,
    chunk_cells: int = 4096,
) -> dict[str, Any]:
    """Sum raw 10x counts into eligible patient × cell-type × site groups.

    The implementation streams columns from the source CSC matrix.  It never
    normalizes before aggregation and reconciles included input counts against
    the output sum.  ``h5py``, NumPy, and SciPy are runtime-only cluster
    dependencies so the repository's metadata-only local environment remains
    lightweight.
    """

    if chunk_cells < 1:
        raise GSE180661Error("chunk_cells must be positive")
    matrix_gate = inspect_10x_h5_matrix(
        matrix_path,
        cell_map_path,
        expected_bytes=expected_bytes,
        frozen_sha256=frozen_sha256,
    )
    try:
        import h5py  # type: ignore[import-not-found]
        import numpy as np  # type: ignore[import-not-found]
        from scipy import sparse  # type: ignore[import-not-found]
    except ImportError as exc:
        raise GSE180661Error(
            "h5py, NumPy, and SciPy are required for pseudobulk aggregation"
        ) from exc

    with group_table_path.open("r", encoding="utf-8", newline="") as handle:
        group_rows = list(csv.DictReader(handle, delimiter="\t"))
    eligible = [row for row in group_rows if row["passes_min_cells"] == "true"]
    if not eligible:
        raise GSE180661Error("no pseudobulk groups pass the pre-specified cell threshold")
    eligible.sort(key=lambda row: row["pseudobulk_group_id"])
    group_ids = [row["pseudobulk_group_id"] for row in eligible]
    group_index = {group_id: index for index, group_id in enumerate(group_ids)}

    cell_to_group: dict[str, int | None] = {}
    with _open_text(cell_map_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            group_id = row["pseudobulk_group_id"]
            cell_to_group[row["cell_id"]] = group_index.get(group_id)

    with h5py.File(matrix_path, "r") as handle:
        matrix = handle["matrix"]
        barcodes = _decode_h5_strings(matrix["barcodes"][:])
        cell_groups = [cell_to_group[barcode] for barcode in barcodes]
        shape = tuple(int(value) for value in matrix["shape"][:])
        n_genes, n_cells = shape
        accumulator = np.zeros((n_genes, len(group_ids)), dtype=np.int64)
        indptr = matrix["indptr"][:]
        included_input_sum = 0
        included_cells = 0
        for cell_start in range(0, n_cells, chunk_cells):
            cell_end = min(n_cells, cell_start + chunk_cells)
            value_start = int(indptr[cell_start])
            value_end = int(indptr[cell_end])
            local_indptr = indptr[cell_start : cell_end + 1] - value_start
            data = matrix["data"][value_start:value_end]
            indices = matrix["indices"][value_start:value_end]
            block = sparse.csc_matrix(
                (data, indices, local_indptr),
                shape=(n_genes, cell_end - cell_start),
                dtype=np.int64,
            )
            selected_local: list[int] = []
            selected_groups: list[int] = []
            for local_index, destination in enumerate(
                cell_groups[cell_start:cell_end]
            ):
                if destination is not None:
                    selected_local.append(local_index)
                    selected_groups.append(destination)
            if not selected_local:
                continue
            assignment = sparse.csr_matrix(
                (
                    np.ones(len(selected_local), dtype=np.int64),
                    (selected_local, selected_groups),
                ),
                shape=(cell_end - cell_start, len(group_ids)),
            )
            partial = (block @ assignment).tocoo()
            np.add.at(
                accumulator,
                (partial.row, partial.col),
                partial.data.astype(np.int64, copy=False),
            )
            included_input_sum += int(block[:, selected_local].sum())
            included_cells += len(selected_local)

        features = matrix["features"]
        gene_ids = _decode_h5_strings(features["id"][:])
        gene_names = _decode_h5_strings(features["name"][:])
        if len(gene_ids) != n_genes or len(gene_names) != n_genes:
            raise GSE180661Error("10x feature vectors do not match matrix gene dimension")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with deterministic_gzip_text_writer(output_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_symbol", *group_ids])
        for index, (gene_id, gene_name) in enumerate(
            zip(gene_ids, gene_names, strict=True)
        ):
            writer.writerow([gene_id, gene_name, *accumulator[index].tolist()])
    output_sum = int(accumulator.sum())
    if output_sum != included_input_sum:
        output_path.unlink(missing_ok=True)
        raise GSE180661Error(
            f"pseudobulk count conservation failed: {included_input_sum} != {output_sum}"
        )
    return {
        "schema_version": "gse180661_pseudobulk_receipt.v1",
        "status": "verified",
        "matrix_gate": matrix_gate,
        "genes": n_genes,
        "matrix_cells": n_cells,
        "included_cells": included_cells,
        "eligible_groups": len(group_ids),
        "included_input_count_sum": included_input_sum,
        "output_count_sum": output_sum,
        "output_path": str(output_path),
        "output_sha256": sha256(output_path),
        "aggregation": "sum_raw_integer_umi_counts",
        "statistical_unit": "patient",
    }
