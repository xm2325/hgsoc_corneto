import csv
import gzip
from pathlib import Path

import pytest

from hgsoc_corneto.external.gse189955 import (
    aggregate_pseudobulk,
    build_corneto_manifest,
    build_group_metadata,
    load_cell_metadata,
    parse_patient_table_rows,
)


def _patient_rows():
    return [
        ["title"],
        [
            "Patient ID",
            "Age at Diagnosis",
            "Disease",
            "Stage",
            "scCOOL-seq2",
            None,
            None,
            None,
            "scRNA-seq",
            None,
            "Number of sampling regions",
            "Anatomical samples",
        ],
        [None] * 12,
        ["OC01", 66, "HGSC", "IIIC", 1, 1, "/", "/", 1, 1, 1, "Primary tumor(1)"],
        [
            "OC21",
            62,
            "Right mesosalpinx cyst & Uterine fibroids",
            "/",
            1,
            1,
            1,
            1,
            "/",
            "/",
            1,
            "Fallopian tube(1)",
        ],
        ["Sum"],
    ]


def _write_metadata(path: Path):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "orig.ident", "nCount_RNA", "nFeature_RNA", "location", "Cell_Type"])
        writer.writerow(["cell_a", "OC01", 10, 4, "PT", "Epithelial cells"])
        writer.writerow(["cell_b", "OC01", 12, 5, "PT", "Epithelial cells"])
        writer.writerow(["cell_c", "OC01", 9, 3, "PT", "Fibroblasts"])
        writer.writerow(["cell_d", "OC21", 8, 3, "FT", "Secretory cells"])


def _write_counts(path: Path, header=None):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header or ["", "cell_a", "cell_b", "cell_c", "cell_d"])
        writer.writerow(["GENE1", 1, 2, 3, 4])
        writer.writerow(["GENE2", 0, 5, 0, 6])


def test_patient_metadata_roles_and_pseudobulk(tmp_path):
    patients = parse_patient_table_rows(_patient_rows())
    assert len(patients) == 2
    assert patients[0]["is_hgsoc"] is True
    assert patients[1]["is_fallopian_tube_reference_donor"] is True

    metadata = tmp_path / "metadata.csv.gz"
    _write_metadata(metadata)
    cells = load_cell_metadata(metadata, patients)
    assert [cell["comparison_role"] for cell in cells] == [
        "hgsoc_epithelial_candidate",
        "hgsoc_epithelial_candidate",
        "hgsoc_site_fibroblast_reference",
        "normal_ft_epithelial_reference",
    ]
    groups = build_group_metadata(cells)
    assert len(groups) == 3
    assert sorted(group["n_cells"] for group in groups) == [1, 1, 2]
    manifest = build_corneto_manifest(groups)
    assert {row["study_accession"] for row in manifest} == {"GSE189955"}
    assert {row["run_accession"] for row in manifest} == {
        row["pseudobulk_group_id"] for row in groups
    }
    assert {row["definitive_malignant"] for row in manifest} == {"false"}
    epithelial = next(
        row for row in manifest if row["comparison_role"] == "hgsoc_epithelial_candidate"
    )
    assert epithelial["malignancy_status"] == "candidate"

    counts = tmp_path / "counts.csv.gz"
    output = tmp_path / "pseudobulk.tsv.gz"
    _write_counts(counts)
    receipt = aggregate_pseudobulk(counts, cells, output, expected_genes=2)
    assert receipt["input_count_sum"] == receipt["output_count_sum"] == 21
    with gzip.open(output, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0][0] == "gene"
    gene1 = dict(zip(rows[0], rows[1], strict=True))
    assert gene1["GSE189955__OC01__primary_tumour__epithelial_cells"] == "3"


def test_count_header_may_reorder_but_must_match_metadata_set(tmp_path):
    patients = parse_patient_table_rows(_patient_rows())
    metadata = tmp_path / "metadata.csv.gz"
    _write_metadata(metadata)
    cells = load_cell_metadata(metadata, patients)
    counts = tmp_path / "counts.csv.gz"
    _write_counts(counts, ["", "cell_b", "cell_a", "cell_c", "cell_d"])
    receipt = aggregate_pseudobulk(counts, cells, tmp_path / "reordered.csv.gz")
    assert receipt["cells"] == 4

    _write_counts(counts, ["", "unknown", "cell_a", "cell_c", "cell_d"])
    with pytest.raises(ValueError, match="do not exactly match"):
        aggregate_pseudobulk(counts, cells, tmp_path / "out.csv.gz")


def test_unknown_location_fails_closed(tmp_path):
    metadata = tmp_path / "metadata.csv.gz"
    with gzip.open(metadata, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "orig.ident", "nCount_RNA", "nFeature_RNA", "location", "Cell_Type"])
        writer.writerow(["cell_x", "OC01", 10, 4, "UNKNOWN", "Epithelial cells"])
    with pytest.raises(ValueError, match="Unknown location code"):
        load_cell_metadata(metadata, parse_patient_table_rows(_patient_rows()))


def test_only_documented_cross_file_cell_id_spellings_are_reconciled(tmp_path):
    metadata = tmp_path / "metadata.csv.gz"
    with gzip.open(metadata, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "orig.ident", "nCount_RNA", "nFeature_RNA", "location", "Cell_Type"])
        writer.writerow(
            ["OC01_PT_Epcam_n_Smart_flow_153", "OC01", 10, 4, "PT", "Epithelial cells"]
        )
        writer.writerow(
            ["OC01_PT_Epcam_p_Smart_mouth_n2", "OC01", 12, 5, "PT", "Fibroblasts"]
        )
    cells = load_cell_metadata(metadata, parse_patient_table_rows(_patient_rows()))
    assert [cell["count_matrix_cell_id"] for cell in cells] == [
        "OC01_PT_Epcam-_Smart_flow_153",
        "OC01_PT_Epcam_p_Smart_mouth-2",
    ]
    assert {cell["cell_id_reconciliation"] for cell in cells} == {
        "documented_geo_file_spelling_difference"
    }
