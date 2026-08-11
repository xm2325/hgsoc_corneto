"""Minimal, deterministic XLSX ingestion for simple source-data worksheets."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from .identifiers import canonicalize_tighe_ocm_id

_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _column_index(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref).group(0)
    value = 0
    for letter in letters:
        value = value * 26 + ord(letter) - ord("A") + 1
    return value - 1


def read_first_worksheet(path: str | Path) -> list[list[str | float | None]]:
    """Read values from the first worksheet using only the Python standard library."""

    with zipfile.ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", _NS):
                shared.append("".join(item.itertext()))

        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        output: list[list[str | float | None]] = []
        for row_node in sheet.findall(".//m:sheetData/m:row", _NS):
            values: dict[int, str | float | None] = {}
            for cell in row_node.findall("m:c", _NS):
                ref = cell.attrib["r"]
                data_type = cell.attrib.get("t")
                value_node = cell.find("m:v", _NS)
                inline_node = cell.find("m:is", _NS)
                raw = value_node.text if value_node is not None else None
                if data_type == "s" and raw is not None:
                    value: str | float | None = shared[int(raw)]
                elif data_type == "inlineStr" and inline_node is not None:
                    value = "".join(inline_node.itertext())
                elif raw is None:
                    value = None
                else:
                    try:
                        value = float(raw)
                    except ValueError:
                        value = raw
                values[_column_index(ref)] = value
            width = max(values, default=-1) + 1
            output.append([values.get(index) for index in range(width)])
        return output


def parse_tighe_abcb1(path: str | Path) -> list[dict[str, object]]:
    """Parse Tighe Table S2 and preserve explicit zero versus missing values."""

    rows = read_first_worksheet(path)
    header_index = next(
        index
        for index, row in enumerate(rows)
        if len(row) >= 2 and str(row[0]).strip() == "OCM" and str(row[1]).strip() == "ABCB1"
    )
    parsed: list[dict[str, object]] = []
    for row in rows[header_index + 1 :]:
        # The published workbook ends with the one-cell footnote
        # ``ND=not done``; it is not an OCM observation.
        if len(row) < 2 or row[0] in (None, "") or row[1] in (None, ""):
            continue
        raw_id = str(row[0]).strip()
        raw_value = row[1] if len(row) > 1 else None
        numeric: float | None
        missing_reason: str | None
        if isinstance(raw_value, (float, int)):
            numeric = float(raw_value)
            missing_reason = None
        else:
            text = "" if raw_value is None else str(raw_value).strip()
            try:
                numeric = float(text)
                missing_reason = None
            except ValueError:
                numeric = None
                missing_reason = text or "blank"
        parsed.append(
            {
                "ocm_id_reported": raw_id,
                "canonical_ocm_id": canonicalize_tighe_ocm_id(raw_id),
                "abcb1_normalized_read_count": numeric,
                "abcb1_missing_reason": missing_reason,
                "source": "Tighe2025 Table S2",
                "source_pmc": "PMC12208324",
            }
        )
    if len(parsed) != 83:
        raise ValueError(f"Expected 83 Table S2 rows, found {len(parsed)}")
    if len({row["canonical_ocm_id"] for row in parsed}) != 83:
        raise ValueError("Duplicate canonical OCM IDs in Table S2")
    return parsed
