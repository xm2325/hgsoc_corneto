from pathlib import Path

import pytest

from hgsoc_corneto.identifiers import read_aliases
from hgsoc_corneto.magetab import load_rna_runs
from hgsoc_corneto.tighe import parse_tighe_table_s1
from hgsoc_corneto.xlsx import parse_tighe_abcb1

ROOT = Path(__file__).resolve().parents[1]


def test_ena_and_sdrf_crosswalk_counts():
    rows = load_rna_runs(
        ROOT / "data/raw/metadata",
        read_aliases(ROOT / "config/ocm_aliases.tsv"),
    )
    assert len(rows) == 117
    assert sum(len(row["fastq_ftp"].split(";")) for row in rows) == 234
    assert sum(row["fastq_total_bytes"] for row in rows) == 477_762_645_114
    assert len({row["run_accession"] for row in rows}) == 117


def test_tighe_abcb1_source_table():
    rows = parse_tighe_abcb1(ROOT / "data/raw/metadata/mmc2.xlsx")
    assert len(rows) == 83
    indexed = {row["canonical_ocm_id"]: row for row in rows}
    assert indexed["OCM328-3"]["abcb1_normalized_read_count"] == 0.0
    assert indexed["OCM326"]["abcb1_missing_reason"] == "ND"


def test_tighe_table_s1_when_supplement_is_cached():
    pdf = ROOT / "tmp/pdfs/mmc1.pdf"
    if not pdf.exists():
        pytest.skip("Run scripts/fetch_public_metadata.py to cache Document S1")
    rows = parse_tighe_table_s1(pdf)
    assert len(rows) == 83
    assert len({row["patient_id"] for row in rows}) == 68
    indexed = {row["canonical_ocm_id"]: row for row in rows}
    assert indexed["OCM231-1"]["chemo_naive_at_biopsy"] is True
    assert indexed["OCM231-5"]["chemo_naive_at_biopsy"] is False
    assert indexed["OCM361a"]["histotype_group"] == "non-HGSOC_or_ambiguous"
