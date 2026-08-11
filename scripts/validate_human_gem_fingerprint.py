#!/usr/bin/env python3
"""Create a read-only structural/provenance receipt for a Human-GEM SBML file."""

from __future__ import annotations

import argparse
import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path

from hgsoc_corneto.metabolic.fingerprint import fingerprint_sbml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sbml", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--expression", type=Path, help="Optional TSV expression matrix for row/column counts")
    parser.add_argument("--manifest", type=Path, help="Optional run manifest for unique study/run counts")
    args = parser.parse_args()
    if not args.sbml.is_file():
        raise FileNotFoundError(args.sbml)
    receipt = fingerprint_sbml(args.sbml)
    metabolite_count = 0
    biomass_ids: list[str] = []
    for _, element in ET.iterparse(args.sbml, events=("end",)):
        tag = element.tag.rsplit("}", 1)[-1]
        if tag == "species":
            metabolite_count += 1
        elif tag == "reaction":
            reaction_id = element.attrib.get("id", "")
            if reaction_id in {"biomass_human", "R_biomass_human"}:
                biomass_ids.append(reaction_id)
        element.clear()
    receipt["metabolites"] = metabolite_count
    receipt["biomass_human_present"] = bool(biomass_ids)
    optional: dict[str, object] = {}
    if args.expression:
        with args.expression.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        optional["expression_rows"] = max(0, len(rows) - 1)
        optional["expression_columns"] = len(rows[0]) if rows else 0
    if args.manifest:
        with args.manifest.open(newline="", encoding="utf-8") as handle:
            manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
        optional["manifest_rows"] = len(manifest_rows)
        optional["manifest_unique_runs"] = len({(r.get("study_accession"), r.get("run_accession")) for r in manifest_rows})
    if optional:
        receipt["optional_input_counts"] = optional
    receipt["biomass_human_present"] = bool(biomass_ids)
    receipt["validation"] = {
        "status": "complete",
        "biomass_detected": bool(receipt["biomass_human_present"]),
        "gpr_partition_complete": bool(receipt["gpr_partition_complete"]),
    }
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
