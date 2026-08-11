#!/usr/bin/env python3
"""Audit pinned Meeson sources without vendoring their unlicensed code."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any

from hgsoc_corneto.metabolic.fingerprint import compare_fingerprint, fingerprint_sbml

EXPECTED_MODEL = {
    "reactions": 13096,
    "genes": 3628,
    "gpr_and": 653,
    "gpr_and_or": 129,
    "gpr_or": 3972,
    "gpr_single_gene": 3282,
    "gpr_no_gene": 5060,
}

KEY_FILES = {
    "integration_implementation": "src/integrate_omics.py",
    "ocm_notebook": "ocm_patient_modelling/Building_OCM-specific_models.ipynb",
    "ocm_flux_output": "ocm_patient_modelling/OCM_clusters_and_fluxes.csv",
    "tpi1_validation": "siTPI1_validation/siTPI1_experimental_simulation_comparison.csv",
}


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, root: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(root)),
        "exists": path.is_file(),
        "bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def notebook_sources(path: Path) -> str:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "code"
    )


def audit_notebook(path: Path, checkout: Path) -> dict[str, Any]:
    source = notebook_sources(path)
    restored = sorted(set(re.findall(r"%store\s+-r\s+([A-Za-z_][A-Za-z0-9_]*)", source)))
    input_csv_references = sorted(
        set(re.findall(r"read_csv\s*\(\s*r?['\"]([^'\"\n]+\.csv)['\"]", source))
    )
    output_csv_references = sorted(
        set(re.findall(r"to_csv\s*\(\s*r?['\"]([^'\"\n]+\.csv)['\"]", source))
    )
    local_absolute = [value for value in input_csv_references if Path(value).is_absolute()]
    unresolved = []
    for value in input_csv_references:
        candidate = Path(value)
        if candidate.is_absolute():
            if not candidate.exists():
                unresolved.append(value)
        elif not (path.parent / candidate).exists() and not (checkout / candidate).exists():
            unresolved.append(value)
    return {
        "restored_ipython_variables": restored,
        "input_csv_references": input_csv_references,
        "output_csv_references": output_csv_references,
        "absolute_input_csv_references": local_absolute,
        "unresolved_input_csv_references": sorted(unresolved),
        "mentions_manual_excel_adjustment": "divide all values by 5 in excel"
        in source.casefold(),
        "mentions_and_or_group_but_no_integration_loop": (
            "ANDORs" in source and "for r in ANDORs" not in source
        ),
    }


def audit_flux_output(path: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = list(reader)
    sample_column = header[0] or "sample"
    samples = [row[0] for row in rows]
    index_rows = [
        {"sample": row[0], "cluster": row[1] if len(row) > 1 else ""} for row in rows
    ]
    reaction_columns = header[2:]
    return (
        {
            "rows": len(rows),
            "columns": len(header),
            "sample_column_as_written": sample_column,
            "cluster_column": header[1] if len(header) > 1 else None,
            "reaction_flux_columns": len(reaction_columns),
            "unique_samples": len(set(samples)),
            "duplicate_samples": sorted(
                sample for sample in set(samples) if samples.count(sample) > 1
            ),
            "contains_biomass_human": "biomass_human" in reaction_columns,
        },
        index_rows,
    )


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "cluster"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--human-gem", type=Path, required=True)
    parser.add_argument("--expected-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    checkout = args.checkout.resolve()
    files = {
        name: file_record(checkout / relative, checkout)
        for name, relative in KEY_FILES.items()
    }
    missing = [name for name, record in files.items() if not record["exists"]]
    if missing:
        raise FileNotFoundError(f"Missing pinned upstream files: {', '.join(missing)}")

    notebook_path = checkout / KEY_FILES["ocm_notebook"]
    flux_path = checkout / KEY_FILES["ocm_flux_output"]
    flux_summary, index_rows = audit_flux_output(flux_path)
    model = fingerprint_sbml(args.human_gem)
    mismatches = compare_fingerprint(model, EXPECTED_MODEL)

    license_candidates = sorted(
        str(path.relative_to(checkout))
        for pattern in ("LICENSE*", "COPYING*")
        for path in checkout.glob(pattern)
        if path.is_file()
    )
    audit = {
        "upstream_repository": "https://github.com/katemeeson/PhD_2024",
        "expected_commit": args.expected_commit,
        "checkout_directory_name": checkout.name,
        "license_files_at_repository_root": license_candidates,
        "vendoring_permitted_by_audit": bool(license_candidates),
        "key_files": files,
        "notebook": audit_notebook(notebook_path, checkout),
        "ocm_flux_output": flux_summary,
        "human_gem": model,
        "human_gem_expected": EXPECTED_MODEL,
        "human_gem_mismatches": mismatches,
        "human_gem_fingerprint_matches": not mismatches,
        "reproduction_boundary": (
            "Exact 49-OCM regeneration is blocked by the absent expression-plus-growth input; "
            "algorithm and public output auditing remain reproducible."
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "upstream_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_tsv(args.output_dir / "ocm_output_index.tsv", index_rows)
    if mismatches:
        raise SystemExit(f"Human-GEM fingerprint mismatch: {mismatches}")


if __name__ == "__main__":
    main()
