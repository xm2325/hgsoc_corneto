import csv
import gzip
import json
from pathlib import Path

import pytest

from hgsoc_corneto.external.gse180661 import (
    CELL_METADATA_COLUMNS,
    GSE180661Error,
    inspect_10x_h5_matrix,
    iter_annotated_cells,
    load_source_manifest,
    parse_mutational_signature_rows,
    parse_sample_inventory_rows,
    write_frozen_matrix_manifest,
)

ROOT = Path(__file__).resolve().parents[1]


def _inventory_rows():
    return [
        [
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
        ],
        [
            "SPECTRUM-OV-003_S1_CD45N_LEFT_ADNEXA",
            "SPECTRUM-OV-003_S1_LEFT_ADNEXA",
            "SPECTRUM-OV-003",
            "hgsoc",
            "Tumor",
            "Adnexa",
            "Left Adnexa",
            "Left Adnexa",
            "Primary",
            "pre-Rx",
            "S1",
            "10x 3' GE",
            "singlet, live, CD45-",
        ],
        [
            "SPECTRUM-OV-003_S1_CD45P_LEFT_ADNEXA",
            "SPECTRUM-OV-003_S1_LEFT_ADNEXA",
            "SPECTRUM-OV-003",
            "hgsoc",
            "Tumor",
            "Adnexa",
            "Left Adnexa",
            "Left Adnexa",
            "Primary",
            "pre-Rx",
            "S1",
            "10x 3' GE",
            "singlet, live, CD45+",
        ],
    ]


def _signature_rows():
    return [
        [
            "patient_id",
            "consensus_signature",
            "consensus_hr_status",
            "wgs_signature",
            "myriad_signature",
            "myriad_gis_score",
        ],
        ["SPECTRUM-OV-003", "HRD-Dup", "HRD", "HRD-Dup", None, None],
    ]


def _cell_row(cell_id: str, sample: str, sort_parameters: str):
    return {
        "cell_id": cell_id,
        "sample": sample,
        "cell_type": "Ovarian.cancer.cell",
        "percent.mt": "5.0",
        "nCount_RNA": "100",
        "nFeature_RNA": "50",
        "umap50_1": "0",
        "umap50_2": "1",
        "cluster_label": "Cancer.1",
        "cluster_label_sub": "NA",
        "cell_type_super": "Ovarian.cancer.super",
        "patient_id": "SPECTRUM-OV-003",
        "tumor_subsite": "Left Adnexa",
        "tumor_site": "Left Adnexa",
        "tumor_supersite": "Adnexa",
        "sort_parameters": sort_parameters,
        "therapy": "pre-Rx",
        "surgery": "S1",
    }


def _write_cells(path: Path, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CELL_METADATA_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def test_manifest_preserves_unfrozen_large_matrix_boundary():
    manifest = load_source_manifest(
        ROOT / "config/external/gse180661_sources.json"
    )
    matrix = next(item for item in manifest["files"] if item["role"] == "count_matrix")
    assert matrix["bytes"] == 32_555_276_423
    assert matrix["sha256"] is None
    assert matrix["default_download"] is False
    assert manifest["expected"]["processed_cells"] == 929_690


def test_cd45_aliquots_collapse_to_same_biological_group(tmp_path):
    inventory = parse_sample_inventory_rows(_inventory_rows())
    signatures = parse_mutational_signature_rows(_signature_rows())
    assert signatures[0]["normalized_signature"] == "hrd_dup"

    metadata = tmp_path / "cells.tsv.gz"
    rows = [
        _cell_row(
            "cell_n",
            "SPECTRUM-OV-003_S1_CD45N_LEFT_ADNEXA",
            "singlet, live, CD45-",
        ),
        _cell_row(
            "cell_p",
            "SPECTRUM-OV-003_S1_CD45P_LEFT_ADNEXA",
            "singlet, live, CD45+",
        ),
    ]
    _write_cells(metadata, rows)
    cells = list(iter_annotated_cells(metadata, inventory, signatures))
    assert {row["pseudobulk_group_id"] for row in cells} == {
        "GSE180661__SPECTRUM-OV-003__ovarian_cancer_cell__left_adnexa"
    }
    assert {row["sort_parameters_reported"] for row in cells} == {
        "singlet, live, CD45-",
        "singlet, live, CD45+",
    }


def test_site_contradiction_fails_closed(tmp_path):
    inventory = parse_sample_inventory_rows(_inventory_rows())
    signatures = parse_mutational_signature_rows(_signature_rows())
    metadata = tmp_path / "cells.tsv.gz"
    row = _cell_row(
        "cell_x",
        "SPECTRUM-OV-003_S1_CD45N_LEFT_ADNEXA",
        "singlet, live, CD45-",
    )
    row["tumor_site"] = "Bowel"
    _write_cells(metadata, [row])
    with pytest.raises(GSE180661Error, match="contradicts Supplementary Table 2"):
        list(iter_annotated_cells(metadata, inventory, signatures))


def test_matrix_gate_requires_digest_frozen_before_schema_read(tmp_path):
    matrix = tmp_path / "matrix.h5"
    matrix.write_bytes(b"\x89HDF\r\n\x1a\n")
    cell_map = tmp_path / "cells.tsv"
    cell_map.write_text("cell_id\ncell_a\n", encoding="utf-8")
    with pytest.raises(GSE180661Error, match="not frozen"):
        inspect_10x_h5_matrix(
            matrix,
            cell_map,
            expected_bytes=8,
            frozen_sha256=None,
        )


def test_unknown_signature_fails_closed():
    rows = _signature_rows()
    rows[1][1] = "HRD-other"
    with pytest.raises(GSE180661Error, match="unknown consensus"):
        parse_mutational_signature_rows(rows)


def test_frozen_matrix_manifest_is_written_as_new_revision(tmp_path):
    source = ROOT / "config/external/gse180661_sources.json"
    output = tmp_path / "frozen.json"
    digest = "a" * 64
    write_frozen_matrix_manifest(
        source,
        output,
        {
            "status": "verified",
            "include_matrix": True,
            "manifest_sha256": "b" * 64,
            "created_at": "2026-08-12T00:00:00+00:00",
            "files": [{"role": "count_matrix", "sha256": digest}],
        },
    )
    frozen = json.loads(output.read_text(encoding="utf-8"))
    matrix = next(item for item in frozen["files"] if item["role"] == "count_matrix")
    assert matrix["sha256"] == digest
    original = load_source_manifest(source)
    original_matrix = next(
        item for item in original["files"] if item["role"] == "count_matrix"
    )
    assert original_matrix["sha256"] is None

    with pytest.raises(GSE180661Error, match="overwrite"):
        write_frozen_matrix_manifest(
            source,
            source,
            {
                "status": "verified",
                "include_matrix": True,
                "files": [{"role": "count_matrix", "sha256": digest}],
            },
        )
