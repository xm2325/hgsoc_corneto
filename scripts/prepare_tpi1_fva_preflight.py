#!/usr/bin/env python3
"""Prepare a solver-free TPI1 deletion/FVA input receipt.

This gate parses Human-GEM SBML directly, validates completed metabolic
baseline receipts, maps the HGNC symbol TPI1 to its model gene product, and
builds the reaction set for a later WT versus gene-deletion FVA.  It does not
import COBRApy/CORNETO/CVXPY and never constructs or solves an optimisation
problem.
"""

from __future__ import annotations

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Any


class PreflightError(ValueError):
    """Raised when an input cannot safely support the planned analysis."""


def _local(value: str) -> str:
    return value.rsplit("}", 1)[-1]


def _attribute(element: ET.Element, local_name: str) -> str | None:
    for key, value in element.attrib.items():
        if _local(key) == local_name:
            return value
    return None


def _reaction_id(sbml_id: str) -> str:
    """Match the reaction ID unescaping used by COBRApy/libSBML."""
    return sbml_id[2:] if sbml_id.startswith("R_") else sbml_id


def _object(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise PreflightError(f"{label} must be an object")
    return value


def _string_list(value: Any, label: str) -> list[str]:
    if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
        raise PreflightError(f"{label} must be an array of non-empty strings")
    return list(value)


def _parse_receipt_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise PreflightError(f"--receipt must be STUDY=PATH, got {spec!r}")
    study, raw_path = spec.split("=", 1)
    if not study or not raw_path:
        raise PreflightError(f"--receipt must be STUDY=PATH, got {spec!r}")
    return study, Path(raw_path)


def _read_receipt(study: str, path: Path) -> dict[str, Any]:
    try:
        root = _object(json.loads(path.read_text(encoding="utf-8")), f"{study} receipt")
    except (OSError, json.JSONDecodeError) as error:
        raise PreflightError(f"{study}: cannot read receipt {path}: {error}") from error
    if root.get("status") != "completed":
        raise PreflightError(f"{study}: baseline status is not completed")
    if root.get("study_accession") != study:
        raise PreflightError(f"{study}: receipt study_accession does not match")

    candidate = _object(root.get("candidate_selection"), f"{study}.candidate_selection")
    candidates = _string_list(candidate.get("selected_reaction_ids"), f"{study} candidates")
    if candidate.get("selected_count") != len(set(candidates)) or len(candidates) != len(set(candidates)):
        raise PreflightError(f"{study}: candidate count or uniqueness check failed")

    corneto = _object(root.get("corneto"), f"{study}.corneto")
    independent = _string_list(corneto.get("independent_active_union"), f"{study} independent union")
    joint = _string_list(corneto.get("joint_active_union"), f"{study} joint union")
    return {
        "path": str(path),
        "sample_count": root.get("sample_count"),
        "candidates": sorted(candidates),
        "independent_active_union": sorted(independent),
        "joint_active_union": sorted(joint),
    }


def _gene_refs(element: ET.Element) -> set[str]:
    refs: set[str] = set()
    for node in element.iter():
        if _local(node.tag) == "geneProductRef":
            gene_id = _attribute(node, "geneProduct")
            if gene_id:
                refs.add(gene_id)
    return refs


def _evaluate_gpr(element: ET.Element, deleted: set[str]) -> bool:
    """Evaluate an SBML FBC association with all non-deleted genes present."""
    tag = _local(element.tag)
    children = list(element)
    if tag == "geneProductRef":
        gene_id = _attribute(element, "geneProduct")
        if not gene_id:
            raise PreflightError("geneProductRef is missing fbc:geneProduct")
        return gene_id not in deleted
    if tag == "and":
        return bool(children) and all(_evaluate_gpr(child, deleted) for child in children)
    if tag == "or":
        return any(_evaluate_gpr(child, deleted) for child in children)
    if tag == "geneProductAssociation":
        if len(children) != 1:
            raise PreflightError("geneProductAssociation must contain exactly one rule")
        return _evaluate_gpr(children[0], deleted)
    raise PreflightError(f"unsupported GPR element {tag!r}")


def _symbol_resources(element: ET.Element) -> set[str]:
    symbols: set[str] = set()
    marker = "identifiers.org/hgnc.symbol/"
    for node in element.iter():
        for value in node.attrib.values():
            if marker in value:
                symbols.add(value.split(marker, 1)[1].split("?", 1)[0].split("#", 1)[0])
    return symbols


def _parse_sbml(path: Path, symbol: str) -> dict[str, Any]:
    try:
        tree = ET.parse(path)
    except (OSError, ET.ParseError) as error:
        raise PreflightError(f"cannot parse Human-GEM SBML: {error}") from error
    root = tree.getroot()

    reaction_elements: dict[str, tuple[str, ET.Element]] = {}
    gene_products: dict[str, dict[str, Any]] = {}
    for element in root.iter():
        tag = _local(element.tag)
        if tag == "geneProduct":
            gene_id = _attribute(element, "id")
            if gene_id:
                gene_products[gene_id] = {
                    "label": _attribute(element, "label"),
                    "name": _attribute(element, "name"),
                    "hgnc_symbols": sorted(_symbol_resources(element)),
                }
        elif tag == "reaction":
            sbml_reaction_id = _attribute(element, "id")
            if sbml_reaction_id:
                reaction_id = _reaction_id(sbml_reaction_id)
                if reaction_id in reaction_elements:
                    raise PreflightError(f"reaction ID normalisation collision for {reaction_id}")
                reaction_elements[reaction_id] = (sbml_reaction_id, element)

    target_gene_ids = sorted(
        gene_id
        for gene_id, metadata in gene_products.items()
        if symbol.casefold() in {value.casefold() for value in metadata["hgnc_symbols"]}
        or symbol.casefold() in {
            str(metadata.get("label") or "").casefold(),
            str(metadata.get("name") or "").casefold(),
            gene_id.casefold(),
        }
    )
    if not target_gene_ids:
        raise PreflightError(f"HGNC symbol {symbol!r} has no Human-GEM gene product mapping")

    all_reaction_ids = set(reaction_elements)
    associated: list[dict[str, Any]] = []
    deleted = set(target_gene_ids)
    for reaction_id, (sbml_reaction_id, reaction) in reaction_elements.items():
        associations = [node for node in reaction if _local(node.tag) == "geneProductAssociation"]
        refs = set().union(*(_gene_refs(node) for node in associations)) if associations else set()
        if not refs.intersection(deleted):
            continue
        before = all(_evaluate_gpr(node, set()) for node in associations)
        after = all(_evaluate_gpr(node, deleted) for node in associations)
        associated.append(
            {
                "reaction_id": reaction_id,
                "sbml_reaction_id": sbml_reaction_id,
                "reaction_name": _attribute(reaction, "name"),
                "gpr_gene_products": sorted(refs),
                "functional_before_deletion": before,
                "functional_after_deletion": after,
                "disabled_by_deletion": before and not after,
                "lower_flux_bound_parameter": _attribute(reaction, "lowerFluxBound"),
                "upper_flux_bound_parameter": _attribute(reaction, "upperFluxBound"),
            }
        )
    if not associated:
        raise PreflightError(f"mapped {symbol} gene product has no associated reactions")
    return {
        "gene_product_count": len(gene_products),
        "reaction_count": len(reaction_elements),
        "all_reaction_ids": all_reaction_ids,
        "target_gene_products": {
            gene_id: gene_products[gene_id] for gene_id in target_gene_ids
        },
        "associated_reactions": sorted(associated, key=lambda row: row["reaction_id"]),
    }


def prepare(
    sbml: Path,
    symbol: str,
    biomass_reaction: str,
    receipt_specs: list[str],
) -> dict[str, Any]:
    if not receipt_specs:
        raise PreflightError("at least one completed baseline --receipt is required")
    parsed_specs = [_parse_receipt_spec(spec) for spec in receipt_specs]
    studies = [study for study, _ in parsed_specs]
    if len(studies) != len(set(studies)):
        raise PreflightError("duplicate study in --receipt")

    model = _parse_sbml(sbml, symbol)
    if biomass_reaction not in model["all_reaction_ids"]:
        raise PreflightError(f"biomass reaction {biomass_reaction!r} is absent from Human-GEM")
    receipts = {study: _read_receipt(study, path) for study, path in parsed_specs}

    provenance: dict[str, set[str]] = defaultdict(set)
    for study, receipt in receipts.items():
        for reaction_id in receipt["candidates"]:
            provenance[reaction_id].add(f"{study}:candidate")
        for reaction_id in receipt["independent_active_union"]:
            provenance[reaction_id].add(f"{study}:independent_active")
        for reaction_id in receipt["joint_active_union"]:
            provenance[reaction_id].add(f"{study}:joint_active")
    for row in model["associated_reactions"]:
        provenance[row["reaction_id"]].add(f"{symbol}:associated")
    provenance[biomass_reaction].add("biomass_objective")

    absent = sorted(set(provenance) - model["all_reaction_ids"])
    if absent:
        raise PreflightError(
            f"{len(absent)} planned FVA reactions are absent from Human-GEM: {absent[:10]}"
        )
    fva_targets = sorted(provenance)
    disabled = sorted(
        row["reaction_id"]
        for row in model["associated_reactions"]
        if row["disabled_by_deletion"]
    )
    if not disabled:
        raise PreflightError(f"deleting {symbol} does not disable any associated reaction")

    return {
        "status": "valid",
        "response_blind": True,
        "solver_called": False,
        "model": {
            "path": str(sbml),
            "reaction_count": model["reaction_count"],
            "gene_product_count": model["gene_product_count"],
            "biomass_reaction": biomass_reaction,
        },
        "gene_deletion": {
            "symbol": symbol,
            "gene_products": model["target_gene_products"],
            "associated_reactions": model["associated_reactions"],
            "disabled_reaction_ids": disabled,
        },
        "baseline_receipts": receipts,
        "planned_fva": {
            "conditions": ["wild_type", f"delta_{symbol}"],
            "reaction_count": len(fva_targets),
            "reaction_ids": fva_targets,
            "reaction_provenance": {
                reaction_id: sorted(values) for reaction_id, values in sorted(provenance.items())
            },
            "not_performed": True,
        },
        "claim_limit": (
            "Input and GPR audit only. WT/deletion growth, FVA ranges, dependency, "
            "and experimental concordance require a later solver-backed analysis."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sbml", type=Path, required=True)
    parser.add_argument("--symbol", default="TPI1")
    parser.add_argument("--biomass-reaction", default="biomass_human")
    parser.add_argument("--receipt", action="append", required=True, metavar="STUDY=PATH")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = prepare(args.sbml, args.symbol, args.biomass_reaction, args.receipt)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    except PreflightError as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": result["status"],
                "fva_reaction_count": result["planned_fva"]["reaction_count"],
                "output": str(args.output),
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
