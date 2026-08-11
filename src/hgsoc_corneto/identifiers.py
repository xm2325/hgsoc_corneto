"""OCM identifier parsing without phenotype-informed guessing."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

CONTROL_NAMES = {"FNE", "KURAMOCHI"}


@dataclass(frozen=True)
class ParsedSourceName:
    raw: str
    source_biospecimen_id: str | None
    patient_id: str | None
    material: str | None
    passage: str | None
    is_control: bool


def _clean(value: str | None) -> str:
    return "" if value is None else " ".join(str(value).strip().split())


def patient_id_from_ocm(ocm_id: str | None) -> str | None:
    """Return the patient-level OCM identifier (for example, ``OCM231``)."""

    if not ocm_id:
        return None
    match = re.match(r"^OCM(\d+)", ocm_id, flags=re.IGNORECASE)
    return f"OCM{match.group(1)}" if match else None


def parse_source_name(raw: str) -> ParsedSourceName:
    """Parse a MAGE-TAB source name into source-level biological identifiers.

    This function standardizes punctuation and separates material/passage. It
    deliberately does not infer omitted biopsy numbers; those mappings live in
    the auditable alias table.
    """

    value = _clean(raw)
    if value.upper() in CONTROL_NAMES:
        return ParsedSourceName(value, None, None, None, None, True)

    if not re.match(r"^OCM(?:[._\-\s]|\d)", value, flags=re.IGNORECASE):
        return ParsedSourceName(value, None, None, None, None, False)

    tail = re.sub(r"^OCM", "", value, flags=re.IGNORECASE)
    tail = tail.lstrip("._- ")
    tail = tail.replace("_", "-").replace(" ", "")
    tail = re.sub(r"\.{1,}", "-", tail)
    tail = re.sub(r"-{2,}", "-", tail).strip("-")

    passage = None
    passage_match = re.search(r"-P(\d+)$", tail, flags=re.IGNORECASE)
    if passage_match:
        passage = f"P{int(passage_match.group(1))}"
        tail = tail[: passage_match.start()]

    tail = re.sub(r"EPCAMNEG", "EPCAMNEG", tail, flags=re.IGNORECASE)
    tail = re.sub(r"EPCAMPOS", "EPCAMPOS", tail, flags=re.IGNORECASE)
    tail = re.sub(r"EPNEG", "EPCAMNEG", tail, flags=re.IGNORECASE)
    tail = re.sub(r"EPPOS", "EPCAMPOS", tail, flags=re.IGNORECASE)

    material = None
    tokens = [token for token in tail.split("-") if token]
    kept: list[str] = []
    for token in tokens:
        upper = token.upper()
        if upper in {"T", "S"} and material is None:
            material = "tumour" if upper == "T" else "stroma"
            continue
        attached = re.fullmatch(r"(\d+)([TS])", upper)
        if attached and material is None:
            kept.append(str(int(attached.group(1))))
            material = "tumour" if attached.group(2) == "T" else "stroma"
            continue
        kept.append(token)

    # Names such as OCM.46-3T-P4 are handled above. A terminal material letter
    # can also be attached to an alphabetic specimen identifier.
    if kept and material is None:
        attached = re.fullmatch(r"(.+?)([TS])", kept[-1], flags=re.IGNORECASE)
        if attached and re.search(r"\d", attached.group(1)):
            kept[-1] = attached.group(1)
            material = "tumour" if attached.group(2).upper() == "T" else "stroma"

    normalized_tail = "-".join(kept).strip("-")
    source_id = f"OCM{normalized_tail}" if normalized_tail else None
    return ParsedSourceName(
        raw=value,
        source_biospecimen_id=source_id,
        patient_id=patient_id_from_ocm(source_id),
        material=material,
        passage=passage,
        is_control=False,
    )


def canonicalize_tighe_ocm_id(raw: str) -> str:
    """Normalize the OCM notation used across Tighe Tables S1 and S2."""

    value = _clean(raw)
    value = re.sub(r"^OCM[.\s_-]*", "", value, flags=re.IGNORECASE)
    value = value.replace(" ", "")
    explicit = {
        "64-3-": "64-3-Ep-",
        "64-3+": "64-3-Ep+",
        "64-3-Ep-": "64-3-Ep-",
        "64-3-Ep+": "64-3-Ep+",
        "376Ta": "376a",
    }
    value = explicit.get(value, value)
    return f"OCM{value}"


def read_aliases(path: str | Path) -> dict[str, str]:
    """Load explicit source-to-canonical mappings and reject duplicates."""

    aliases: dict[str, str] = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source = row["source_biospecimen_id"].strip()
            canonical = row["canonical_ocm_id"].strip()
            if source in aliases and aliases[source] != canonical:
                raise ValueError(f"Conflicting alias for {source}")
            aliases[source] = canonical
    return aliases


def apply_alias(source_biospecimen_id: str | None, aliases: dict[str, str]) -> str | None:
    if source_biospecimen_id is None:
        return None
    return aliases.get(source_biospecimen_id, source_biospecimen_id)
