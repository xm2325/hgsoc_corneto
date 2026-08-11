#!/usr/bin/env python3
"""Fail-closed, solver-free gate for a future TPI1/FVA analysis.

The script checks Human-GEM XML gene/reaction identifiers and, when present,
selected reaction IDs from CORNETO receipts.  It does not load COBRApy,
CORNETO, CVXPY, or any solver license, and it never performs a deletion or
flux-range calculation.
"""

from __future__ import annotations

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1].casefold()


def _gene_tokens(element: ET.Element) -> set[str]:
    tokens: set[str] = set()
    for key, value in element.attrib.items():
        if _local(key) in {"id", "name", "label"} and value:
            tokens.add(value.strip())
    for child in element.iter():
        if child is element:
            continue
        if _local(child.tag) in {"geneproductref", "gene", "geneproduct"}:
            for key, value in child.attrib.items():
                if _local(key) in {"id", "geneproduct", "geneproductid", "name", "label"} and value:
                    tokens.add(value.strip())
    return tokens


def _reaction_id(sbml_id: str) -> str:
    """Match the reaction ID unescaping used by COBRApy/libSBML."""
    return sbml_id[2:] if sbml_id.startswith("R_") else sbml_id


def _symbol_resources(element: ET.Element) -> set[str]:
    symbols: set[str] = set()
    marker = "identifiers.org/hgnc.symbol/"
    for node in element.iter():
        for value in node.attrib.values():
            if marker in value:
                symbols.add(value.split(marker, 1)[1].split("?", 1)[0].split("#", 1)[0])
    return symbols


def _parse_model(path: Path) -> tuple[set[str], set[str], list[str]]:
    genes: set[str] = set()
    reactions: set[str] = set()
    target_gene_ids: set[str] = set()
    reaction_elements: list[ET.Element] = []
    try:
        root = ET.parse(path).getroot()
        for element in root.iter():
            tag = _local(element.tag)
            attrs = {(_local(key), value.strip()) for key, value in element.attrib.items() if value}
            if tag in {"geneproduct", "gene"}:
                ids = {value for key, value in attrs if key in {"id", "name", "label"}}
                genes.update(ids)
                if any("tpi1" in value.casefold() for value in ids | _symbol_resources(element)):
                    target_gene_ids.update(ids)
            elif tag == "reaction":
                ids = {value for key, value in attrs if key == "id"}
                reactions.update(_reaction_id(value) for value in ids)
                reaction_elements.append(element)
    except (OSError, ET.ParseError) as error:
        raise RuntimeError(f"cannot parse SBML: {error}") from error
    tpi1_reactions: set[str] = set()
    for element in reaction_elements:
        if _gene_tokens(element).intersection(target_gene_ids):
            reaction_id = next(
                (
                    _reaction_id(value)
                    for key, value in element.attrib.items()
                    if _local(key) == "id"
                ),
                None,
            )
            if reaction_id:
                tpi1_reactions.add(reaction_id)
    tpi1_hits = [f"gene:{value}" for value in sorted(target_gene_ids)]
    tpi1_hits.extend(f"reaction:{value}" for value in sorted(tpi1_reactions))
    return genes, reactions, tpi1_hits


def _receipt_reactions(paths: list[Path]) -> tuple[set[str], list[str]]:
    selected: set[str] = set()
    missing: list[str] = []
    for path in paths:
        try:
            root: Any = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            missing.append(str(path))
            continue
        candidate = root.get("candidate_selection") if isinstance(root, dict) else None
        ids = candidate.get("selected_reaction_ids") if isinstance(candidate, dict) else None
        if not isinstance(ids, list):
            missing.append(str(path))
            continue
        selected.update(value for value in ids if isinstance(value, str) and value)
    return selected, missing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sbml", type=Path, required=True)
    parser.add_argument("--receipt", action="append", type=Path, default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        genes, reactions, tpi1_hits = _parse_model(args.sbml)
        selected, missing_receipts = _receipt_reactions(args.receipt)
        absent_selected = sorted(selected - reactions)
        result = {
            "status": "valid" if tpi1_hits and not absent_selected and not missing_receipts else "blocked",
            "response_blind": True,
            "sbml": str(args.sbml),
            "reaction_count": len(reactions),
            "gene_count": len(genes),
            "tpi1_hits": tpi1_hits,
            "selected_reaction_count": len(selected),
            "selected_reactions_absent_from_model": absent_selected,
            "missing_or_invalid_receipts": missing_receipts,
            "next_step": "Run WT/Delta-TPI1 deletion and FVA only after this gate is valid; this script performs neither.",
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.output.with_suffix(args.output.suffix + ".tmp")
        tmp.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        tmp.replace(args.output)
    except RuntimeError as error:
        print(str(error), file=sys.stderr)
        return 1
    print(json.dumps({"status": result["status"], "output": str(args.output)}))
    return 0 if result["status"] == "valid" else 2


if __name__ == "__main__":
    sys.exit(main())
