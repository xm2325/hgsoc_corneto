"""Dependency-light fingerprints for SBML models.

The parser deliberately uses the standard library so that provenance checks can
run before COBRApy or CORNETO is installed. It counts FBC gene-association tags
inside each reaction and therefore does not depend on how a particular SBML
reader normalizes identifiers.
"""

from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _gpr_class(reaction: ET.Element) -> str:
    tags = [_local_name(element.tag) for element in reaction.iter()]
    gene_count = tags.count("geneProductRef")
    has_and = "and" in tags
    has_or = "or" in tags
    if gene_count == 0:
        return "gpr_no_gene"
    if gene_count == 1:
        return "gpr_single_gene"
    if has_and and has_or:
        return "gpr_and_or"
    if has_and:
        return "gpr_and"
    if has_or:
        return "gpr_or"
    # An association with multiple refs but no logical operator is malformed.
    return "gpr_unclassified"


def _sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_sbml(path: str | Path) -> dict[str, Any]:
    """Return a stable structural fingerprint for an SBML model."""

    model_path = Path(path)
    counts = {
        "reactions": 0,
        "genes": 0,
        "gpr_and": 0,
        "gpr_and_or": 0,
        "gpr_or": 0,
        "gpr_single_gene": 0,
        "gpr_no_gene": 0,
        "gpr_unclassified": 0,
    }
    biomass_sbml_ids: list[str] = []

    for _, element in ET.iterparse(model_path, events=("end",)):
        name = _local_name(element.tag)
        if name == "geneProduct":
            counts["genes"] += 1
        elif name == "reaction":
            counts["reactions"] += 1
            counts[_gpr_class(element)] += 1
            reaction_id = element.attrib.get("id", "")
            if "biomass" in reaction_id.casefold():
                biomass_sbml_ids.append(reaction_id)
            element.clear()

    classified = sum(value for key, value in counts.items() if key.startswith("gpr_"))
    return {
        "path": str(model_path),
        "bytes": model_path.stat().st_size,
        "sha256": _sha256(model_path),
        **counts,
        "gpr_partition_complete": classified == counts["reactions"],
        "biomass_sbml_ids": biomass_sbml_ids,
    }


def compare_fingerprint(
    observed: dict[str, Any], expected: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    """Return expected/observed mismatches for keys present in ``expected``."""

    return {
        key: {"expected": expected_value, "observed": observed.get(key)}
        for key, expected_value in expected.items()
        if observed.get(key) != expected_value
    }
