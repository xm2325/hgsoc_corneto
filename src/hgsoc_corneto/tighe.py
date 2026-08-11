"""Extraction of Tighe et al. 2025 supplementary patient/OCM metadata."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .identifiers import canonicalize_tighe_ocm_id


def _text(value: Any) -> str | None:
    if value in (None, ""):
        return None
    return " ".join(str(value).replace("\n", " ").split())


def _coalesce(row: list[Any], *indices: int) -> str | None:
    for index in indices:
        if index < len(row):
            value = _text(row[index])
            if value is not None:
                return value
    return None


def _histotype_group(reported: str | None) -> str | None:
    if reported is None:
        return None
    return "HGSOC" if reported.strip() == "HGSOC" else "non-HGSOC_or_ambiguous"


def _parse_page_rows(page_number: int, table: list[list[Any]]) -> list[dict[str, Any]]:
    parsed: list[dict[str, Any]] = []
    current_patient: dict[str, Any] = {}
    previous_ocm: dict[str, Any] = {}

    for row in table:
        if page_number == 22:
            if not (_coalesce(row, 0) or _coalesce(row, 10)):
                continue
            number = _coalesce(row, 0)
            if number and not number.isdigit():
                continue
            patient_values = {
                "table_row": number,
                "patient_numeric": _coalesce(row, 2),
                "figo_stage": _coalesce(row, 3),
                "histotype_reported": _coalesce(row, 4),
                "primary_tp53": _coalesce(row, 5, 7),
            }
            ocm_values = {
                "ocm_id_reported": _coalesce(row, 10),
                "chemo_naive_at_biopsy": _coalesce(row, 12),
                "biopsy_type": _coalesce(row, 13),
                "ocm_tp53_dna": _coalesce(row, 14),
                "ocm_tp53_protein": _coalesce(row, 16),
                "p53_without_nutlin3": _coalesce(row, 18),
                "p53_with_nutlin3": _coalesce(row, 19),
                "references": _coalesce(row, 21),
            }
        elif page_number == 23:
            patient_values = {
                "table_row": _coalesce(row, 0),
                "patient_numeric": _coalesce(row, 1),
                "figo_stage": _coalesce(row, 2),
                "histotype_reported": _coalesce(row, 3),
                "primary_tp53": _coalesce(row, 4, 5),
            }
            ocm_values = {
                "ocm_id_reported": _coalesce(row, 7),
                "chemo_naive_at_biopsy": _coalesce(row, 8),
                "biopsy_type": _coalesce(row, 9),
                "ocm_tp53_dna": _coalesce(row, 10),
                "ocm_tp53_protein": _coalesce(row, 11),
                "p53_without_nutlin3": _coalesce(row, 12),
                "p53_with_nutlin3": _coalesce(row, 13),
                "references": _coalesce(row, 14),
            }
        elif page_number == 24:
            patient_values = {
                "table_row": _coalesce(row, 0),
                "patient_numeric": _coalesce(row, 1),
                "figo_stage": _coalesce(row, 2),
                "histotype_reported": _coalesce(row, 3),
                "primary_tp53": _coalesce(row, 4),
            }
            ocm_values = {
                "ocm_id_reported": _coalesce(row, 5),
                "chemo_naive_at_biopsy": _coalesce(row, 6),
                "biopsy_type": _coalesce(row, 7),
                "ocm_tp53_dna": _coalesce(row, 8),
                "ocm_tp53_protein": _coalesce(row, 9),
                "p53_without_nutlin3": _coalesce(row, 10),
                "p53_with_nutlin3": _coalesce(row, 11),
                "references": _coalesce(row, 12),
            }
        else:
            raise ValueError(f"Unexpected page {page_number}")

        if patient_values["patient_numeric"]:
            current_patient = patient_values
        if not ocm_values["ocm_id_reported"]:
            continue
        if not current_patient:
            raise ValueError(f"OCM row without patient context on page {page_number}")

        # In the published table, the chemo-naive and biopsy-type cells for the
        # OCM64-3 EpCAM fractions are vertically merged. Fill only these two
        # design fields on continuation rows; never carry molecular values.
        if not patient_values["patient_numeric"] and previous_ocm:
            for field in ("chemo_naive_at_biopsy", "biopsy_type"):
                if ocm_values[field] is None:
                    ocm_values[field] = previous_ocm.get(field)

        record = {**current_patient, **ocm_values}
        record["patient_id"] = f"OCM{record['patient_numeric']}"
        record["canonical_ocm_id"] = canonicalize_tighe_ocm_id(record["ocm_id_reported"])
        record["histotype_group"] = _histotype_group(record["histotype_reported"])
        record["chemo_naive_at_biopsy"] = {
            "Y": True,
            "N": False,
        }.get(record["chemo_naive_at_biopsy"])
        record["source"] = "Tighe2025 Table S1"
        record["source_pmc"] = "PMC12208324"
        record["source_pdf_page"] = page_number
        parsed.append(record)
        previous_ocm = ocm_values
    return parsed


def parse_tighe_table_s1(pdf_path: str | Path) -> list[dict[str, Any]]:
    """Extract Table S1 from pages 22–24 of Document S1.

    The parser is intentionally pinned to the published table layout and fails
    loudly if row or patient counts change.
    """

    try:
        import pdfplumber
    except ImportError as exc:  # pragma: no cover - dependency error path
        raise RuntimeError("Install the metadata extra to parse the Tighe PDF") from exc

    records: list[dict[str, Any]] = []
    with pdfplumber.open(pdf_path) as pdf:
        if len(pdf.pages) < 24:
            raise ValueError("Tighe Document S1 must contain at least 24 pages")
        page22_tables = pdf.pages[21].extract_tables()
        page23_tables = pdf.pages[22].extract_tables()
        page24_tables = pdf.pages[23].extract_tables()
        if not page22_tables or not page23_tables or not page24_tables:
            raise ValueError("Could not locate Table S1 ruling lines")
        records.extend(_parse_page_rows(22, page22_tables[0]))
        records.extend(_parse_page_rows(23, page23_tables[0]))
        records.extend(_parse_page_rows(24, page24_tables[0]))

    if len(records) != 83:
        raise ValueError(f"Expected 83 OCM rows, found {len(records)}")
    if len({record["patient_id"] for record in records}) != 68:
        raise ValueError("Expected 68 unique patients")
    if len({record["canonical_ocm_id"] for record in records}) != 83:
        raise ValueError("Duplicate canonical OCM IDs in Table S1")
    return records
