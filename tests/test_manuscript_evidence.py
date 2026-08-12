from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts import build_manuscript_evidence as evidence

ROOT = Path(__file__).resolve().parents[1]


def _manifest_rows() -> list[dict[str, str]]:
    with (ROOT / "data/processed/metadata/ocm_master_manifest.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_frozen_manifest_statistics() -> None:
    stats = evidence.derive_manifest_statistics(_manifest_rows())
    assert stats["manifest_rows"] == 117
    assert stats["primary_ocms"] == 60
    assert stats["primary_patients"] == 52
    assert stats["patient_multiplicity"] == {"1": 45, "2": 6, "3": 1}
    assert stats["repeated_patients"] == 7
    assert stats["rows_from_repeated_patients"] == 15
    assert stats["one_per_patient_selection_count"] == 192
    assert [stats["by_study"][study]["primary_hgsoc_tumour"] for study in evidence.STUDY_ORDER] == [
        9,
        13,
        11,
        27,
    ]
    assert stats["fastq_objects"] == 234
    assert stats["fastq_bytes"] == 477_762_645_114


def test_dispositions_are_exhaustive() -> None:
    counts: dict[str, int] = {}
    for row in _manifest_rows():
        disposition = evidence.classify_manifest_row(row)
        counts[disposition] = counts.get(disposition, 0) + 1
    assert counts == {
        "primary_hgsoc_tumour": 60,
        "stroma_reference": 33,
        "cell_line_control": 2,
        "tumour_non_hgsoc_or_ambiguous": 17,
        "tumour_not_in_tighe_screen": 4,
        "tumour_nonrepresentative_duplicate": 1,
    }


def test_latex_escape() -> None:
    assert evidence.latex_escape("OCM64_3 & 5%") == r"OCM64\_3 \& 5\%"


def test_claim_status_gate_rejects_unknown() -> None:
    row = {
        "claim_id": "X",
        "question": "q",
        "estimand": "e",
        "falsification_rule": "f",
        "status": "proven",
        "evidence_source": "s",
        "permitted_wording": "p",
        "prohibited_wording": "n",
    }
    with pytest.raises(ValueError, match="invalid claim status"):
        evidence._validate_claims([row])
