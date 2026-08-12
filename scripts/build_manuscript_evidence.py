#!/usr/bin/env python3
"""Build deterministic manuscript evidence tables from frozen inputs.

The script never contacts the network and never reads solver licences. Remote
model values enter only through ``evidence/roihu_result_snapshot.json``, which
records the authoritative receipt path, SHA256 and extraction boundary.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from typing import Any

ALLOWED_CLAIM_STATUSES = {
    "supported under tested conditions",
    "weakened",
    "falsified",
    "blocked",
    "pending",
}
STUDY_ORDER = ("E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _truth(value: str) -> bool:
    return value.strip().lower() == "true"


def _display(value: str) -> str:
    return "" if value in {"NA", "None", "null"} else value


def classify_manifest_row(row: dict[str, str]) -> str:
    """Return the mutually exclusive frozen analysis disposition."""
    if _truth(row["primary_cohort_eligible"]):
        return "primary_hgsoc_tumour"
    if row["sample_class"] == "stroma":
        return "stroma_reference"
    if row["sample_class"] == "cell_line_control":
        return "cell_line_control"
    if row["sample_class"] == "tumour":
        if not _truth(row["in_tighe_83_ocm_screen"]):
            return "tumour_not_in_tighe_screen"
        if row["histotype_group"] != "HGSOC":
            return "tumour_non_hgsoc_or_ambiguous"
        if not _truth(row["is_representative_rna_library"]):
            return "tumour_nonrepresentative_duplicate"
        return "tumour_excluded_other"
    return "unclassified"


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def _atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False, newline=""
    ) as handle:
        handle.write(content)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    _atomic_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _atomic_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False, newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def _validate_claims(rows: list[dict[str, str]]) -> None:
    seen: set[str] = set()
    for row in rows:
        claim_id = row.get("claim_id", "")
        if not claim_id or claim_id in seen:
            raise ValueError(f"invalid or duplicate claim_id: {claim_id!r}")
        seen.add(claim_id)
        if row.get("status") not in ALLOWED_CLAIM_STATUSES:
            raise ValueError(f"invalid claim status for {claim_id}: {row.get('status')!r}")
        for key in (
            "question",
            "estimand",
            "falsification_rule",
            "evidence_source",
            "permitted_wording",
            "prohibited_wording",
        ):
            if not row.get(key):
                raise ValueError(f"missing {key} for {claim_id}")


def _validate_remote_snapshot(payload: dict[str, Any]) -> None:
    if payload.get("schema_version") != "manuscript_roihu_receipt_snapshot.v1":
        raise ValueError("unexpected Roihu snapshot schema")
    receipts = payload.get("receipts")
    if not isinstance(receipts, dict) or not receipts:
        raise ValueError("Roihu snapshot has no receipts")
    for name, receipt in receipts.items():
        if not isinstance(receipt, dict):
            raise ValueError(f"receipt {name} is not an object")
        if name != "rna_aggregation":
            digest = receipt.get("sha256")
            if not isinstance(digest, str) or len(digest) != 64:
                raise ValueError(f"receipt {name} has invalid SHA256")
            if not receipt.get("path"):
                raise ValueError(f"receipt {name} has no path")


def derive_manifest_statistics(rows: list[dict[str, str]]) -> dict[str, Any]:
    if len(rows) != 117:
        raise ValueError(f"expected 117 manifest rows, observed {len(rows)}")
    run_ids = [row["run_accession"] for row in rows]
    if len(set(run_ids)) != len(run_ids):
        raise ValueError("run_accession is not unique")
    if set(row["study_accession"] for row in rows) != set(STUDY_ORDER):
        raise ValueError("unexpected study set")

    dispositions = Counter(classify_manifest_row(row) for row in rows)
    primary = [row for row in rows if classify_manifest_row(row) == "primary_hgsoc_tumour"]
    patient_counts = Counter(row["patient_id"] for row in primary)
    multiplicity = Counter(patient_counts.values())
    fastq_objects = 0
    fastq_bytes = 0
    by_study: dict[str, dict[str, Any]] = {}
    for study in STUDY_ORDER:
        subset = [row for row in rows if row["study_accession"] == study]
        study_primary = [
            row for row in subset if classify_manifest_row(row) == "primary_hgsoc_tumour"
        ]
        study_dispositions = Counter(classify_manifest_row(row) for row in subset)
        study_fastq = []
        for row in subset:
            study_fastq.extend(
                int(item) for item in row["fastq_bytes"].split(";") if item and item != "NA"
            )
        fastq_objects += len(study_fastq)
        fastq_bytes += sum(study_fastq)
        by_study[study] = {
            "all_runs": len(subset),
            "primary_hgsoc_tumour": len(study_primary),
            "primary_patients_within_study": len({row["patient_id"] for row in study_primary}),
            "stroma_reference": study_dispositions["stroma_reference"],
            "cell_line_control": study_dispositions["cell_line_control"],
            "tumour_non_hgsoc_or_ambiguous": study_dispositions["tumour_non_hgsoc_or_ambiguous"],
            "tumour_not_in_tighe_screen": study_dispositions["tumour_not_in_tighe_screen"],
            "tumour_nonrepresentative_duplicate": study_dispositions[
                "tumour_nonrepresentative_duplicate"
            ],
            "fastq_objects": len(study_fastq),
            "fastq_bytes": sum(study_fastq),
        }

    expected_dispositions = {
        "primary_hgsoc_tumour": 60,
        "stroma_reference": 33,
        "cell_line_control": 2,
        "tumour_non_hgsoc_or_ambiguous": 17,
        "tumour_not_in_tighe_screen": 4,
        "tumour_nonrepresentative_duplicate": 1,
    }
    for name, expected in expected_dispositions.items():
        if dispositions[name] != expected:
            raise ValueError(f"unexpected {name}: {dispositions[name]} != {expected}")
    if dispositions["tumour_excluded_other"] or dispositions["unclassified"]:
        raise ValueError(f"unexpected residual dispositions: {dict(dispositions)}")
    if len(primary) != 60 or len(patient_counts) != 52:
        raise ValueError("primary OCM/patient invariant failed")
    if dict(sorted(multiplicity.items())) != {1: 45, 2: 6, 3: 1}:
        raise ValueError(f"unexpected patient multiplicity: {dict(multiplicity)}")
    if fastq_objects != 234 or fastq_bytes != 477_762_645_114:
        raise ValueError(f"FASTQ invariant failed: {fastq_objects}, {fastq_bytes}")

    return {
        "manifest_rows": len(rows),
        "primary_ocms": len(primary),
        "primary_patients": len(patient_counts),
        "repeated_patients": sum(count > 1 for count in patient_counts.values()),
        "rows_from_repeated_patients": sum(count for count in patient_counts.values() if count > 1),
        "one_per_patient_selection_count": 2**6 * 3,
        "patient_multiplicity": {str(key): value for key, value in sorted(multiplicity.items())},
        "dispositions": dict(sorted(dispositions.items())),
        "fastq_objects": fastq_objects,
        "fastq_bytes": fastq_bytes,
        "fastq_gib": fastq_bytes / (1024**3),
        "by_study": by_study,
    }


def _registry_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    order = {study: index for index, study in enumerate(STUDY_ORDER)}
    for row in sorted(
        rows, key=lambda item: (order[item["study_accession"]], item["run_accession"])
    ):
        output.append(
            {
                "study_accession": row["study_accession"],
                "ena_project": row["ena_project"],
                "run_accession": row["run_accession"],
                "ena_sample_accession": row["ena_sample_accession"],
                "raw_source_name": row["source_name"],
                "canonical_ocm_id": row["canonical_ocm_id"],
                "patient_id": row["patient_id"],
                "sample_class": row["sample_class"],
                "passage": _display(row["passage"]),
                "histotype_reported": _display(row["histotype_reported"]),
                "histotype_group": _display(row["histotype_group"]),
                "tighe_table_row": _display(row["tighe_table_row"]),
                "representative_rna_library": row["is_representative_rna_library"],
                "primary_cohort_eligible": row["primary_cohort_eligible"],
                "disposition": classify_manifest_row(row),
                "first_public": row["first_public"],
                "metadata_source": row["metadata_source"],
            }
        )
    return output


def _numeric_ledger(
    stats: dict[str, Any],
    remote: dict[str, Any],
    tighe_rows: list[dict[str, str]],
    tighe_abcb1_rows: list[dict[str, str]],
    meeson_audit: dict[str, Any],
    toy_global: dict[str, Any],
    toy_joint: dict[str, Any],
) -> list[dict[str, object]]:
    """Return a tidy source/derivation ledger for manuscript result numbers."""
    rows: list[dict[str, object]] = []

    def add(
        evidence_id: str,
        value: object,
        unit: str,
        evidence_class: str,
        source_path: str,
        source_locator: str,
        derivation: str,
        claim_limit: str,
    ) -> None:
        rows.append(
            {
                "evidence_id": evidence_id,
                "value": value,
                "unit": unit,
                "evidence_class": evidence_class,
                "source_path": source_path,
                "source_locator": source_locator,
                "derivation": derivation,
                "claim_limit": claim_limit,
            }
        )

    manifest_source = "data/processed/metadata/ocm_master_manifest.tsv"
    for name, value, unit, derivation in (
        ("cohort.manifest_runs", stats["manifest_rows"], "runs", "count unique run_accession"),
        ("cohort.primary_ocms", stats["primary_ocms"], "OCMs", "count primary_hgsoc_tumour rows"),
        (
            "cohort.primary_patients",
            stats["primary_patients"],
            "patients",
            "count distinct patient_id among primary rows",
        ),
        (
            "cohort.repeated_patients",
            stats["repeated_patients"],
            "patients",
            "count primary patient_id values represented by more than one row",
        ),
        (
            "cohort.rows_from_repeated_patients",
            stats["rows_from_repeated_patients"],
            "OCMs",
            "sum rows contributed by repeated patients",
        ),
        (
            "cohort.one_per_patient_selections",
            stats["one_per_patient_selection_count"],
            "valid selections",
            "product of repeated-family sizes: 2^6 times 3",
        ),
        (
            "cohort.fastq_objects",
            stats["fastq_objects"],
            "files",
            "split each fastq_bytes field on semicolon and count",
        ),
        (
            "cohort.fastq_bytes",
            stats["fastq_bytes"],
            "bytes",
            "sum every semicolon-delimited fastq_bytes integer",
        ),
    ):
        add(
            name,
            value,
            unit,
            "deterministic derivation",
            manifest_source,
            "all rows",
            derivation,
            "Input/cohort accounting only",
        )
    for study in STUDY_ORDER:
        item = stats["by_study"][study]
        for field, unit in (
            ("all_runs", "runs"),
            ("primary_hgsoc_tumour", "OCMs"),
            ("primary_patients_within_study", "patients"),
            ("stroma_reference", "runs"),
            ("cell_line_control", "runs"),
            ("tumour_non_hgsoc_or_ambiguous", "runs"),
            ("tumour_not_in_tighe_screen", "runs"),
            ("tumour_nonrepresentative_duplicate", "runs"),
            ("fastq_objects", "files"),
            ("fastq_bytes", "bytes"),
        ):
            add(
                f"cohort.{study}.{field}",
                item[field],
                unit,
                "deterministic derivation",
                manifest_source,
                f"study_accession={study}",
                f"group by study and compute {field}",
                "Input/cohort accounting only",
            )

    receipts = remote["receipts"]
    for study in receipts["rna_aggregation"]["studies"]:
        for field, unit in (
            ("run_count", "runs"),
            ("gene_count", "genes"),
            ("transcript_count", "transcripts"),
        ):
            add(
                f"rna_aggregation.{study['study_accession']}.{field}",
                study[field],
                unit,
                "scientific receipt",
                study["path"],
                f"/{field}; sha256={study['sha256']}",
                "direct JSON field",
                receipts["rna_aggregation"]["claim_limit"],
            )

    rank_summary = receipts["nmf_pooled_rank_summary"]
    for rank, record in rank_summary["rows"].items():
        for field in (
            "sample_count",
            "gene_count",
            "nmf_runs",
            "cophenetic_correlation",
            "average_silhouette",
            "converged_runs",
        ):
            add(
                f"nmf.pooled.rank{rank}.{field}",
                record[field],
                "count"
                if field.endswith("count") or field in {"nmf_runs", "converged_runs"}
                else "metric",
                "model receipt",
                rank_summary["path"],
                f"rank={rank}; column={field}; sha256={rank_summary['sha256']}",
                "direct TSV field",
                rank_summary["claim_limit"],
            )
    nmf_compare = receipts["nmf_pooled_vs_cohort"]
    for record in nmf_compare["cohort_metrics"]:
        for field in (
            "sample_count",
            "adjusted_rand_index",
            "normalized_mutual_information",
            "mapped_assignment_agreement",
            "mean_matched_loading_spearman",
        ):
            add(
                f"nmf.{record['study_accession']}.{field}",
                record[field],
                "count" if field == "sample_count" else "metric",
                "model receipt",
                nmf_compare["path"],
                (
                    f"/cohort_metrics/{record['study_accession']}/{field}; "
                    f"sha256={nmf_compare['sha256']}"
                ),
                "direct JSON field",
                nmf_compare["claim_limit"],
            )
    for field in ("chi_square", "p_value", "cramers_v"):
        add(
            f"nmf.state_by_study.{field}",
            nmf_compare["pooled_state_by_study"][field],
            "metric",
            "model receipt",
            nmf_compare["path"],
            f"/pooled_state_by_study/{field}; sha256={nmf_compare['sha256']}",
            "direct JSON field",
            nmf_compare["claim_limit"],
        )
    balanced_nmf = receipts["nmf_patient_balanced"]
    for field in (
        "common_run_count",
        "adjusted_rand_index",
        "normalized_mutual_information",
        "mapped_assignment_agreement",
    ):
        add(
            f"nmf.patient_balanced.{field}",
            balanced_nmf[field],
            "count" if field == "common_run_count" else "metric",
            "model receipt",
            balanced_nmf["path"],
            f"/{field}; sha256={balanced_nmf['sha256']}",
            "direct JSON field",
            balanced_nmf["claim_limit"],
        )

    grid = receipts["regulatory_grid"]
    add(
        "regulatory.grid.fit_count",
        grid["fit_count"],
        "fits",
        "model receipt",
        grid["path"],
        f"/fit_count; sha256={grid['sha256']}",
        "five modes times nine lambdas, receipt-validated",
        grid["claim_limit"],
    )
    add(
        "regulatory.grid.lambda_count",
        len(grid["lambda_values"]),
        "values",
        "configured and receipt-validated",
        grid["path"],
        f"/lambda_values; sha256={grid['sha256']}",
        "length of frozen lambda array",
        grid["claim_limit"],
    )
    regulatory_balanced = receipts["regulatory_patient_balanced"]
    for field in (
        "lambda_nominal",
        "common_run_count",
        "pooled_edge_union",
        "patient_balanced_edge_union",
        "union_jaccard",
        "mean_sample_jaccard",
    ):
        add(
            f"regulatory.patient_balanced.{field}",
            regulatory_balanced[field],
            "metric",
            "model receipt",
            regulatory_balanced["path"],
            f"/{field}; sha256={regulatory_balanced['sha256']}",
            "direct JSON field",
            regulatory_balanced["claim_limit"],
        )
    richer = receipts["regulatory_richer_policy"]
    for group in richer["groups"]:
        for field in ("sample_count", "union_jaccard", "mean_sample_jaccard"):
            add(
                f"regulatory.richer.{group['group']}.{field}",
                group[field],
                "metric",
                "model receipt",
                richer["path"],
                f"/groups/{group['group']}/{field}; sha256={richer['sha256']}",
                "direct JSON field",
                richer["claim_limit"],
            )
    alternatives = receipts["regulatory_alternative_optima"]
    for field in (
        "expected_samples",
        "usable_nonempty_samples",
        "median_accepted_alternative_count",
        "median_core_edge_count",
        "median_edge_union_count",
        "median_incumbent_edge_count",
        "median_mean_pairwise_jaccard",
    ):
        add(
            f"regulatory.alternative_optima.{field}",
            alternatives[field],
            "metric",
            "model receipt",
            alternatives["path"],
            f"/{field}; sha256={alternatives['sha256']}",
            "direct JSON field",
            alternatives["claim_limit"],
        )
    fingerprint = receipts["human_gem_fingerprint"]
    for field, unit in (
        ("bytes", "bytes"),
        ("reactions", "reactions"),
        ("metabolites", "metabolites"),
        ("genes", "genes"),
    ):
        add(
            f"human_gem.{field}",
            fingerprint[field],
            unit,
            "model fingerprint receipt",
            fingerprint["path"],
            f"/{field}; sha256={fingerprint['sha256']}",
            "direct JSON field",
            "Model identity only; no OCM-specific result",
        )
    pooled_gate = receipts["pooled_metabolic_input_gate"]
    for field, unit in (
        ("sample_count", "samples"),
        ("ocm_count", "OCMs"),
        ("patient_count", "patients"),
        ("gene_count", "genes"),
    ):
        add(
            f"metabolic.pooled_input.{field}",
            pooled_gate[field],
            unit,
            "input-gate receipt",
            pooled_gate["path"],
            f"/{field}; sha256={pooled_gate['sha256']}",
            "direct JSON field",
            pooled_gate["claim_limit"],
        )
    add(
        "tighe.table_s1.ocms",
        len(tighe_rows),
        "OCMs",
        "frozen supplementary table",
        "data/processed/metadata/tighe_ocm_characteristics.tsv",
        "data rows; hash in manuscript_evidence_snapshot.json",
        "count canonical_ocm_id rows",
        "Wider biobank screen, not the primary RNA set",
    )
    add(
        "tighe.table_s1.patients",
        len({row["patient_id"] for row in tighe_rows}),
        "patients",
        "frozen supplementary table",
        "data/processed/metadata/tighe_ocm_characteristics.tsv",
        "patient_id; one per patient after de-duplication",
        "count distinct patient_id",
        "Wider biobank screen, not the primary RNA set",
    )
    add(
        "tighe.table_s2.ocms",
        len(tighe_abcb1_rows),
        "OCMs",
        "frozen supplementary table",
        "data/processed/metadata/tighe_abcb1.tsv",
        "data rows",
        "count rows",
        "ABCB1 covariate table; not an exact AUC/GI50 matrix",
    )
    flux = meeson_audit["ocm_flux_output"]
    for field, unit in (("rows", "OCM labels"), ("reaction_flux_columns", "reaction-flux columns")):
        add(
            f"meeson.public_output.{field}",
            flux[field],
            unit,
            "public-output audit receipt",
            "data/processed/meeson/upstream_audit.json",
            f"/ocm_flux_output/{field}; hash in manuscript_evidence_snapshot.json",
            "direct JSON field",
            meeson_audit["reproduction_boundary"],
        )
    add(
        "meeson.toy.retained_set_jaccard",
        toy_global["sequential"]["retained_set_jaccard"],
        "Jaccard",
        "toy model receipt",
        "data/processed/meeson/toy_corneto_global_retention.json",
        "/sequential/retained_set_jaccard",
        "direct JSON field",
        "Toy mechanism test only",
    )
    add(
        "meeson.toy.global_growth",
        toy_global["global"]["primary"]["growth"],
        "model growth",
        "toy model receipt",
        "data/processed/meeson/toy_corneto_global_retention.json",
        "/global/primary/growth",
        "direct JSON field",
        "Toy mechanism test only",
    )
    add(
        "meeson.toy.global_objective",
        toy_global["global"]["primary"]["objective_value"],
        "objective",
        "toy model receipt",
        "data/processed/meeson/toy_corneto_global_retention.json",
        "/global/primary/objective_value",
        "direct JSON field",
        "Toy mechanism test only",
    )
    for field in ("independent_active_union_size", "joint_active_union_size"):
        add(
            f"meeson.toy.{field}",
            toy_joint["comparison"][field],
            "reactions",
            "toy model receipt",
            "data/processed/meeson/toy_corneto_joint_fba.json",
            f"/comparison/{field}",
            "direct JSON field",
            "Toy mechanism test only; smaller union is encouraged by the objective",
        )
    return rows


def _cohort_tex(stats: dict[str, Any]) -> str:
    lines = ["% Generated by scripts/build_manuscript_evidence.py; do not edit."]
    for study in STUDY_ORDER:
        item = stats["by_study"][study]
        exclusions = []
        for count, label in (
            (item["stroma_reference"], "stroma"),
            (item["cell_line_control"], "cell-line controls"),
            (item["tumour_non_hgsoc_or_ambiguous"], "non-HGSOC or ambiguous tumour"),
            (item["tumour_not_in_tighe_screen"], "tumour runs not matched to Tighe Table S1"),
            (
                item["tumour_nonrepresentative_duplicate"],
                "non-representative duplicate tumour library",
            ),
        ):
            if count:
                exclusions.append(f"{count} {label}")
        lines.append(
            f"{latex_escape(study)} & {item['all_runs']} & {item['primary_hgsoc_tumour']} & "
            f"{item['primary_patients_within_study']} & {latex_escape('; '.join(exclusions))} \\\\"
        )
    lines.append("\\midrule")
    lines.append("Total & 117 & 60 & 52 (global) & 57 non-primary/reference runs \\\\")
    return "\n".join(lines) + "\n"


def _values_tex(stats: dict[str, Any], remote: dict[str, Any]) -> str:
    nmf = remote["receipts"]["nmf_pooled_rank_summary"]["rows"]
    return "\n".join(
        [
            "% Generated by scripts/build_manuscript_evidence.py; do not edit.",
            f"\\newcommand{{\\ManifestRunCount}}{{{stats['manifest_rows']}}}",
            f"\\newcommand{{\\PrimaryOCMCount}}{{{stats['primary_ocms']}}}",
            f"\\newcommand{{\\PrimaryPatientCount}}{{{stats['primary_patients']}}}",
            f"\\newcommand{{\\FastqObjectCount}}{{{stats['fastq_objects']}}}",
            f"\\newcommand{{\\FastqByteCount}}{{{stats['fastq_bytes']:,}}}",
            f"\\newcommand{{\\FastqGiB}}{{{stats['fastq_gib']:.3f}}}",
            f"\\newcommand{{\\NmfRankTwoCophenetic}}{{{nmf['2']['cophenetic_correlation']:.3f}}}",
            f"\\newcommand{{\\NmfRankTwoSilhouette}}{{{nmf['2']['average_silhouette']:.3f}}}",
            f"\\newcommand{{\\NmfRankThreeCophenetic}}{{{nmf['3']['cophenetic_correlation']:.3f}}}",
            f"\\newcommand{{\\NmfRankThreeSilhouette}}{{{nmf['3']['average_silhouette']:.3f}}}",
            "",
        ]
    )


def _nmf_tex(remote: dict[str, Any]) -> str:
    rows = remote["receipts"]["nmf_pooled_vs_cohort"]["cohort_metrics"]
    lines = ["% Generated by scripts/build_manuscript_evidence.py; do not edit."]
    for row in rows:
        lines.append(
            f"{latex_escape(row['study_accession'])} & {row['adjusted_rand_index']:.3f} & "
            f"{row['normalized_mutual_information']:.3f} & "
            f"{row['mapped_assignment_agreement']:.3f} & "
            f"{row['mean_matched_loading_spearman']:.3f} \\\\"
        )
    balanced = remote["receipts"]["nmf_patient_balanced"]
    lines.append("\\midrule")
    lines.append(
        f"Patient-balanced 52 & {balanced['adjusted_rand_index']:.3f} & "
        f"{balanced['normalized_mutual_information']:.3f} & "
        f"{balanced['mapped_assignment_agreement']:.3f} & -- \\\\"
    )
    return "\n".join(lines) + "\n"


def _claim_tex(rows: list[dict[str, str]]) -> str:
    lines = ["% Generated by scripts/build_manuscript_evidence.py; do not edit."]
    for row in rows:
        lines.append(
            "{} & {} & {} & {} & {} \\\\".format(
                latex_escape(row["claim_id"]),
                latex_escape(row["question"]),
                latex_escape(row["estimand"]),
                latex_escape(row["falsification_rule"]),
                latex_escape(row["status"]),
            )
        )
    return "\n".join(lines) + "\n"


def _failure_tex(rows: list[dict[str, str]]) -> str:
    lines = ["% Generated by scripts/build_manuscript_evidence.py; do not edit."]
    for row in rows:
        lines.append(
            "{} & {} & {} & {} & {} \\\\".format(
                latex_escape(row["failure_id"]),
                latex_escape(row["analysis_class"]),
                latex_escape(row["observed_event"]),
                latex_escape(row["correction_or_disposition"]),
                latex_escape(row["scientific_impact"]),
            )
        )
    return "\n".join(lines) + "\n"


def _paper_tex(rows: list[dict[str, str]]) -> str:
    lines = ["% Generated by scripts/build_manuscript_evidence.py; do not edit."]
    for row in rows:
        identifiers = row["ocm_ids"].replace(";", ", ")
        lines.append(
            "{} ({}) & {} & {} & {} & {} \\\\".format(
                latex_escape(row["paper_key"]),
                latex_escape(row["year"]),
                latex_escape(row["reported_scope"]),
                latex_escape(identifiers),
                latex_escape(row["source_locator"]),
                latex_escape(row["claim_limit"]),
            )
        )
    return "\n".join(lines) + "\n"


def _registry_tex(rows: list[dict[str, object]]) -> str:
    lines = ["% Generated by scripts/build_manuscript_evidence.py; do not edit."]
    for row in rows:
        lines.append(
            "{} & {} & {} & {} & {} & {} & {} & {} \\\\".format(
                latex_escape(row["study_accession"]),
                latex_escape(row["run_accession"]),
                latex_escape(row["raw_source_name"]),
                latex_escape(row["canonical_ocm_id"]),
                latex_escape(row["patient_id"]),
                latex_escape(row["sample_class"]),
                latex_escape(row["histotype_group"]),
                latex_escape(row["disposition"]),
            )
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()

    manifest_path = root / "data/processed/metadata/ocm_master_manifest.tsv"
    manifest_audit_path = root / "data/processed/metadata/manifest_audit.json"
    tighe_path = root / "data/processed/metadata/tighe_ocm_characteristics.tsv"
    tighe_abcb1_path = root / "data/processed/metadata/tighe_abcb1.tsv"
    meeson_audit_path = root / "data/processed/meeson/upstream_audit.json"
    meeson_index_path = root / "data/processed/meeson/ocm_output_index.tsv"
    toy_global_path = root / "data/processed/meeson/toy_corneto_global_retention.json"
    toy_joint_path = root / "data/processed/meeson/toy_corneto_joint_fba.json"
    claims_path = root / "evidence/claims.tsv"
    failures_path = root / "evidence/failures.tsv"
    papers_path = root / "evidence/paper_ocm_evidence.tsv"
    remote_path = root / "evidence/roihu_result_snapshot.json"

    manifest = _read_tsv(manifest_path)
    tighe_rows = _read_tsv(tighe_path)
    tighe_abcb1_rows = _read_tsv(tighe_abcb1_path)
    meeson_index_rows = _read_tsv(meeson_index_path)
    claims = _read_tsv(claims_path)
    failures = _read_tsv(failures_path)
    papers = _read_tsv(papers_path)
    remote = json.loads(remote_path.read_text(encoding="utf-8"))
    meeson_audit = json.loads(meeson_audit_path.read_text(encoding="utf-8"))
    toy_global = json.loads(toy_global_path.read_text(encoding="utf-8"))
    toy_joint = json.loads(toy_joint_path.read_text(encoding="utf-8"))
    _validate_claims(claims)
    _validate_remote_snapshot(remote)
    stats = derive_manifest_statistics(manifest)
    registry = _registry_rows(manifest)

    if len(tighe_rows) != 83:
        raise ValueError("Tighe Table S1 row count is not 83")
    if len(tighe_abcb1_rows) != 83:
        raise ValueError("Tighe Table S2 row count is not 83")
    tighe_paper = next(row for row in papers if row["paper_key"] == "tighe2025")
    if tighe_paper["ocm_ids"].split(";") != [row["canonical_ocm_id"] for row in tighe_rows]:
        raise ValueError("paper evidence Tighe OCM IDs do not exactly match Table S1")
    for paper_key, reference_number in (
        ("pillay2019", "1"),
        ("nelson2020", "2"),
        ("barnes2021", "3"),
        ("coulsongilmer2021", "4"),
    ):
        expected_ids = [
            row["canonical_ocm_id"]
            for row in tighe_rows
            if reference_number in re.findall(r"\d+", row["references"])
        ]
        observed_ids = next(row for row in papers if row["paper_key"] == paper_key)[
            "ocm_ids"
        ].split(";")
        if observed_ids != expected_ids:
            raise ValueError(f"{paper_key} crosswalk does not match Tighe references")
    meeson_resource = next(row for row in papers if row["paper_key"] == "meeson_phd_2024")
    if meeson_resource["ocm_ids"].split(";") != [row["sample"] for row in meeson_index_rows]:
        raise ValueError("paper evidence does not exactly preserve Meeson output labels")

    registry_fields = [
        "study_accession",
        "ena_project",
        "run_accession",
        "ena_sample_accession",
        "raw_source_name",
        "canonical_ocm_id",
        "patient_id",
        "sample_class",
        "passage",
        "histotype_reported",
        "histotype_group",
        "tighe_table_row",
        "representative_rna_library",
        "primary_cohort_eligible",
        "disposition",
        "first_public",
        "metadata_source",
    ]
    _atomic_tsv(root / "evidence/study_ocm_registry.tsv", registry_fields, registry)
    numeric_fields = [
        "evidence_id",
        "value",
        "unit",
        "evidence_class",
        "source_path",
        "source_locator",
        "derivation",
        "claim_limit",
    ]
    _atomic_tsv(
        root / "evidence/numeric_ledger.tsv",
        numeric_fields,
        _numeric_ledger(
            stats, remote, tighe_rows, tighe_abcb1_rows, meeson_audit, toy_global, toy_joint
        ),
    )

    source_paths = [
        manifest_path,
        manifest_audit_path,
        tighe_path,
        tighe_abcb1_path,
        meeson_audit_path,
        meeson_index_path,
        toy_global_path,
        toy_joint_path,
        claims_path,
        failures_path,
        papers_path,
        remote_path,
    ]
    for study in STUDY_ORDER:
        source_paths.extend(
            [
                root / f"data/raw/metadata/biostudies/{study}.sdrf.txt",
                root / f"data/raw/metadata/biostudies/{study}.idf.txt",
            ]
        )
    snapshot = {
        "schema_version": "manuscript_evidence_snapshot.v1",
        "derivation_script": "scripts/build_manuscript_evidence.py",
        "source_sha256": {str(path.relative_to(root)): _sha256(path) for path in source_paths},
        "cohort": stats,
        "claim_status_counts": dict(sorted(Counter(row["status"] for row in claims).items())),
        "claim_count": len(claims),
        "failure_count": len(failures),
        "paper_resource_count": len(papers),
        "remote_snapshot_sha256": _sha256(remote_path),
    }
    _atomic_json(root / "evidence/manuscript_evidence_snapshot.json", snapshot)
    _atomic_text(root / "tex/generated/evidence_values.tex", _values_tex(stats, remote))
    _atomic_text(root / "tex/generated/cohort_table_rows.tex", _cohort_tex(stats))
    _atomic_text(root / "tex/generated/nmf_table_rows.tex", _nmf_tex(remote))
    _atomic_text(root / "tex/generated/claim_ledger_rows.tex", _claim_tex(claims))
    _atomic_text(root / "tex/generated/failure_ledger_rows.tex", _failure_tex(failures))
    _atomic_text(root / "tex/generated/paper_ocm_evidence_rows.tex", _paper_tex(papers))
    _atomic_text(root / "tex/generated/study_ocm_registry_rows.tex", _registry_tex(registry))

    print(
        json.dumps(
            {
                "status": "generated",
                "runs": stats["manifest_rows"],
                "primary_ocms": stats["primary_ocms"],
                "primary_patients": stats["primary_patients"],
                "claims": len(claims),
                "papers_resources": len(papers),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
