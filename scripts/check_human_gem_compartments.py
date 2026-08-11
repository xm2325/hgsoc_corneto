#!/usr/bin/env python3
"""Inventory Human-GEM compartments to replace heuristic external selection."""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1].casefold()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sbml", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    compartments: list[dict[str, object]] = []
    for _, element in ET.iterparse(args.sbml, events=("end",)):
        if _local(element.tag) == "compartment":
            attrs = {key.rsplit("}", 1)[-1].casefold(): value for key, value in element.attrib.items()}
            compartments.append({
                "id": attrs.get("id", ""),
                "name": attrs.get("name", ""),
                "outside": attrs.get("outside", ""),
                "boundary_condition": attrs.get("boundarycondition", ""),
            })
        element.clear()
    explicit = [
        item for item in compartments
        if any(token in f"{item['id']} {item['name']}".casefold() for token in ("extracellular", "external", "outside"))
    ]
    result = {
        "status": "valid" if explicit else "blocked",
        "compartment_count": len(compartments),
        "explicit_external_candidates": explicit,
        "all_compartments": compartments,
        "claim_limit": "inventory only; no flux solve performed",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "compartment_count": len(compartments), "output": str(args.output)}))
    return 0 if explicit else 2


if __name__ == "__main__":
    raise SystemExit(main())
