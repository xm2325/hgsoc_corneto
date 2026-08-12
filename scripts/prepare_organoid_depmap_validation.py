#!/usr/bin/env python3
"""Audit GSE208216 or extract a provenance-scoped DepMap candidate subset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from hgsoc_corneto.external_validation import (
    audit_gse_count_matrix,
    depmap_download_preflight,
    download_verified,
    extract_depmap_gene_effects,
    extract_gse_candidate_log_cpm,
    file_sha256,
    iter_candidate_symbols,
    read_csv_rows,
    read_gene_map,
    read_tsv_rows,
    resolve_depmap_hgsoc_models,
    summarize_dependency_rows,
    summarize_expression_rows,
)
from hgsoc_corneto.io import write_json, write_tsv


def _load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def audit_gse(args: argparse.Namespace) -> None:
    contract = _load_json(args.contract)
    expected_samples = [row["sample_id"] for row in contract["samples"]]
    result = audit_gse_count_matrix(
        args.counts,
        expected_sha256=contract["count_matrix"]["sha256"],
        expected_samples=expected_samples,
        expected_gene_rows=contract["count_matrix"]["gene_rows"],
    )
    result.update(
        {
            "accession": contract["accession"],
            "source_url": contract["count_matrix"]["url"],
            "sample_metadata": contract["samples"],
            "claim_limits": [contract["scope_note"], contract["copy_number_note"]],
        }
    )
    write_json(args.receipt, result)


def prepare_gse(args: argparse.Namespace) -> None:
    contract = _load_json(args.contract)
    samples = contract["samples"]
    audit = audit_gse_count_matrix(
        args.counts,
        expected_sha256=contract["count_matrix"]["sha256"],
        expected_samples=[row["sample_id"] for row in samples],
        expected_gene_rows=contract["count_matrix"]["gene_rows"],
    )
    candidates = sorted(set(iter_candidate_symbols(args.candidates)))
    expression, missing_genes = extract_gse_candidate_log_cpm(
        args.counts,
        sample_groups={row["sample_id"]: row["sample_class"] for row in samples},
        gene_id_to_symbol=read_gene_map(args.gene_map),
        candidate_symbols=candidates,
    )
    if missing_genes:
        raise ValueError(f"candidate genes absent after explicit Ensembl mapping: {missing_genes}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "candidate_expression_log_cpm.tsv", expression)
    write_tsv(
        args.output_dir / "candidate_expression_summary.tsv",
        summarize_expression_rows(expression),
    )
    write_json(
        args.output_dir / "receipt.json",
        {
            **audit,
            "accession": contract["accession"],
            "candidate_gene_count": len(candidates),
            "candidate_value_count": len(expression),
            "missing_candidate_genes": missing_genes,
            "gene_map": {"path": str(args.gene_map), "sha256": file_sha256(args.gene_map)},
            "candidates": {
                "path": str(args.candidates),
                "sha256": file_sha256(args.candidates),
            },
            "claim_limit": (
                "The 11 PDO versus 3 fallopian-tube organoid summary is descriptive external "
                "model evidence; it is not a patient cohort or a clinical contrast."
            ),
        },
    )


def fetch_gse(args: argparse.Namespace) -> None:
    contract = _load_json(args.contract)
    source = contract["count_matrix"]
    target = args.output_dir / Path(source["url"]).name
    download = download_verified(
        source["url"],
        target,
        expected_sha256=source["sha256"],
    )
    audit = audit_gse_count_matrix(
        target,
        expected_sha256=source["sha256"],
        expected_samples=[row["sample_id"] for row in contract["samples"]],
        expected_gene_rows=source["gene_rows"],
    )
    write_json(
        args.receipt,
        {
            **audit,
            "accession": contract["accession"],
            "source_url": source["url"],
            "download": download,
            "sample_metadata": contract["samples"],
            "claim_limits": [contract["scope_note"], contract["copy_number_note"]],
        },
    )


def depmap_preflight(args: argparse.Namespace) -> None:
    contract = _load_json(args.contract)
    receipt = depmap_download_preflight(
        release=args.release,
        model_path=args.models,
        gene_effect_path=args.gene_effect,
        release_readme_path=args.release_readme,
        landing_url=contract["official_download_landing_page"],
    )
    write_json(args.receipt, receipt)


def prepare_depmap(args: argparse.Namespace) -> None:
    model_rows = read_csv_rows(args.models)
    curated_rows = read_tsv_rows(args.curated_models)
    selected, comparators = resolve_depmap_hgsoc_models(model_rows, curated_rows)
    model_groups = {
        row["ModelID"]: row["validation_group"] for row in [*selected, *comparators]
    }
    candidates = sorted(set(iter_candidate_symbols(args.candidates)))
    effects, missing_models, missing_genes = extract_depmap_gene_effects(
        args.gene_effect,
        model_groups=model_groups,
        candidate_symbols=candidates,
    )
    if missing_models and not args.allow_unscreened_models:
        raise ValueError(
            "selected/comparator models missing from gene-effect matrix; rerun with "
            f"--allow-unscreened-models only to preserve them as missing: {missing_models}"
        )
    if missing_genes:
        raise ValueError(f"candidate genes absent from gene-effect matrix: {missing_genes}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "depmap_models.tsv", [*selected, *comparators])
    write_tsv(args.output_dir / "candidate_gene_effects.tsv", effects)
    write_tsv(
        args.output_dir / "candidate_dependency_summary.tsv",
        summarize_dependency_rows(effects),
    )
    write_json(
        args.output_dir / "receipt.json",
        {
            "scientific_success": bool(effects) and not missing_genes,
            "release": args.release,
            "inputs": {
                "Model.csv": {"path": str(args.models), "sha256": file_sha256(args.models)},
                "CRISPRGeneEffect.csv": {
                    "path": str(args.gene_effect),
                    "sha256": file_sha256(args.gene_effect),
                },
                "candidates": {
                    "path": str(args.candidates),
                    "sha256": file_sha256(args.candidates),
                },
            },
            "curated_hgsoc_like_model_ids": [row["ModelID"] for row in selected],
            "other_ovarian_model_ids": [row["ModelID"] for row in comparators],
            "missing_screened_models": missing_models,
            "missing_candidate_genes": missing_genes,
            "candidate_gene_count": len(candidates),
            "extracted_value_count": len(effects),
            "claim_limit": (
                "Chronos gene effect is in-vitro functional plausibility, not clinical drug "
                "response; other ovarian models are comparators, not HGSOC positives."
            ),
        },
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    gse = commands.add_parser("audit-gse208216")
    gse.add_argument("--counts", type=Path, required=True)
    gse.add_argument(
        "--contract",
        type=Path,
        default=Path("config/external_validation/gse208216.json"),
    )
    gse.add_argument("--receipt", type=Path, required=True)
    gse.set_defaults(func=audit_gse)

    gse_fetch = commands.add_parser("fetch-gse208216")
    gse_fetch.add_argument("--output-dir", type=Path, required=True)
    gse_fetch.add_argument("--receipt", type=Path, required=True)
    gse_fetch.add_argument(
        "--contract",
        type=Path,
        default=Path("config/external_validation/gse208216.json"),
    )
    gse_fetch.set_defaults(func=fetch_gse)

    gse_prepare = commands.add_parser("prepare-gse208216")
    gse_prepare.add_argument("--counts", type=Path, required=True)
    gse_prepare.add_argument("--gene-map", type=Path, required=True)
    gse_prepare.add_argument("--candidates", type=Path, required=True)
    gse_prepare.add_argument("--output-dir", type=Path, required=True)
    gse_prepare.add_argument(
        "--contract",
        type=Path,
        default=Path("config/external_validation/gse208216.json"),
    )
    gse_prepare.set_defaults(func=prepare_gse)

    depmap = commands.add_parser("prepare-depmap")
    depmap.add_argument("--models", type=Path, required=True)
    depmap.add_argument("--gene-effect", type=Path, required=True)
    depmap.add_argument("--candidates", type=Path, required=True)
    depmap.add_argument("--release", required=True)
    depmap.add_argument("--output-dir", type=Path, required=True)
    depmap.add_argument(
        "--curated-models",
        type=Path,
        default=Path("config/external_validation/depmap_hgsoc_models.tsv"),
    )
    depmap.add_argument("--allow-unscreened-models", action="store_true")
    depmap.set_defaults(func=prepare_depmap)

    depmap_gate = commands.add_parser("depmap-preflight")
    depmap_gate.add_argument("--release")
    depmap_gate.add_argument("--models", type=Path)
    depmap_gate.add_argument("--gene-effect", type=Path)
    depmap_gate.add_argument("--release-readme", type=Path)
    depmap_gate.add_argument("--receipt", type=Path, required=True)
    depmap_gate.add_argument(
        "--contract",
        type=Path,
        default=Path("config/external_validation/depmap_download_contract.json"),
    )
    depmap_gate.set_defaults(func=depmap_preflight)
    return root


def main() -> None:
    args = parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
