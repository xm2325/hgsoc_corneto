#!/usr/bin/env python3
"""Validate and prepare the label-blind Taylor/HGSOC analysis design.

This script only reads metadata.  It never reads, derives, or fills paclitaxel
phenotype values.  The outputs are design manifests and deterministic
patient-grouped fold assignments for the cross-sectional, intrinsic, and
acquired tracks.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


MISSING = frozenset({"", "NA", "N/A", "null", "None"})
REQUIRED_MANIFEST_FIELDS = frozenset(
    {
        "study_accession",
        "run_accession",
        "canonical_ocm_id",
        "patient_id",
        "sample_class",
        "histotype_group",
        "is_representative_rna_library",
        "hgsoc_tumour_eligible",
        "primary_cohort_eligible",
        "exact_paclitaxel_auc_available",
        "exact_paclitaxel_gi50_available",
        "exact_cumulative_paclitaxel_exposure_available",
        "paclitaxel_auc",
        "paclitaxel_gi50_nm",
        "cumulative_paclitaxel_exposure_mg",
    }
)
REQUIRED_FAMILY_FIELDS = frozenset(
    {
        "patient_id",
        "ocm_ids",
        "n_ocms",
        "relationship_type",
        "analysis_role",
    }
)
OUTCOME_SPECS = {
    "paclitaxel_auc": {
        "available": "exact_paclitaxel_auc_available",
        "kind": "nonnegative",
    },
    "paclitaxel_gi50_nm": {
        "available": "exact_paclitaxel_gi50_available",
        "kind": "positive",
    },
    "cumulative_paclitaxel_exposure_mg": {
        "available": "exact_cumulative_paclitaxel_exposure_available",
        "kind": "nonnegative",
    },
}


def _is_missing(value: str | None) -> bool:
    return value is None or value.strip() in MISSING


def _read_tsv(path: Path, required: set[str] | frozenset[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        fields = list(reader.fieldnames)
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate TSV columns in {path}")
        missing = sorted(required - set(fields))
        if missing:
            raise ValueError(f"{path} is missing required columns: {missing}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"TSV has no data rows: {path}")
    return rows


def _boolean(row: dict[str, str], field: str, *, context: str) -> bool:
    value = row.get(field, "").strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"{context}: {field} must be true/false, got {row.get(field)!r}")


def _number(
    row: dict[str, str],
    field: str,
    *,
    context: str,
    kind: str,
) -> float | None:
    value = row.get(field, "")
    if _is_missing(value):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{context}: {field} is not numeric: {value!r}") from error
    if not math.isfinite(parsed):
        raise ValueError(f"{context}: {field} must be finite: {value!r}")
    if kind == "nonnegative" and parsed < 0:
        raise ValueError(f"{context}: {field} must be non-negative: {parsed}")
    if kind == "positive" and parsed <= 0:
        raise ValueError(f"{context}: {field} must be positive: {parsed}")
    return parsed


def _validate_outcomes(rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    availability: dict[str, dict[str, int]] = {}
    for field, spec in OUTCOME_SPECS.items():
        available_field = str(spec["available"])
        counts = {"available": 0, "missing": 0}
        for row in rows:
            context = f"{row.get('study_accession')}/{row.get('run_accession')}"
            flag = _boolean(row, available_field, context=context)
            value = _number(row, field, context=context, kind=str(spec["kind"]))
            if flag and value is None:
                raise ValueError(
                    f"{context}: {available_field}=true but {field} is missing"
                )
            if not flag and value is not None:
                raise ValueError(
                    f"{context}: {available_field}=false but {field} is populated; "
                    "do not create phenotype values in the metadata manifest"
                )
            counts["available" if flag else "missing"] += 1
        availability[field] = counts
    return availability


def _primary_rows(rows: list[dict[str, str]], schema: dict[str, Any]) -> list[dict[str, str]]:
    cohort = schema["cohort"]
    primary = [
        row
        for row in rows
        if _boolean(row, "primary_cohort_eligible", context=row["run_accession"])
    ]
    expected = int(cohort["primary_runs"])
    if len(primary) != expected:
        raise ValueError(f"Expected {expected} primary rows, found {len(primary)}")
    required_sample_class = str(cohort["required_sample_class"])
    required_histotype = str(cohort["required_histotype_group"])
    required_representative = str(cohort["required_representative_flag"])
    for row in primary:
        context = row["run_accession"]
        if row["sample_class"] != required_sample_class:
            raise ValueError(f"{context}: primary row is not tumour")
        if row["histotype_group"] != required_histotype:
            raise ValueError(f"{context}: primary row is not HGSOC")
        if row["is_representative_rna_library"] != required_representative:
            raise ValueError(f"{context}: primary row is not the representative RNA library")
        if not _boolean(row, "hgsoc_tumour_eligible", context=context):
            raise ValueError(f"{context}: primary row is not hgsoc_tumour_eligible")
        for field in ("canonical_ocm_id", "patient_id"):
            if _is_missing(row.get(field)):
                raise ValueError(f"{context}: primary row has missing {field}")
    if len({row["canonical_ocm_id"] for row in primary}) != int(cohort["primary_ocms"]):
        raise ValueError("Primary cohort OCM count does not match the frozen schema")
    if len({row["patient_id"] for row in primary}) != int(cohort["primary_patients"]):
        raise ValueError("Primary cohort patient count does not match the frozen schema")
    return primary


def _assign_group_folds(
    rows: list[dict[str, str]], *, n_splits: int
) -> tuple[dict[str, int], list[int]]:
    """Assign complete patient groups to balanced folds deterministically."""

    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["patient_id"]].append(row)
    if len(groups) < n_splits:
        raise ValueError(f"Need at least {n_splits} patient groups, found {len(groups)}")

    loads = [0] * n_splits
    assignment: dict[str, int] = {}
    ordered = sorted(groups.items(), key=lambda item: (-len(item[1]), item[0]))
    for patient_id, members in ordered:
        fold = min(range(n_splits), key=lambda index: (loads[index], index))
        assignment[patient_id] = fold
        loads[fold] += len(members)
    return assignment, loads


def _family_design(
    primary: list[dict[str, str]], families: list[dict[str, str]]
) -> tuple[list[dict[str, str]], dict[str, str]]:
    by_ocm = {row["canonical_ocm_id"]: row for row in primary}
    ocm_to_family: dict[str, str] = {}
    design: list[dict[str, str]] = []
    for family in families:
        if family["relationship_type"] != "longitudinal":
            continue
        family_id = family["patient_id"]
        ids = [value for value in family["ocm_ids"].split(";") if value]
        try:
            declared_count = int(family["n_ocms"])
        except ValueError as error:
            raise ValueError(f"Invalid n_ocms for longitudinal family {family_id}") from error
        if declared_count != len(ids):
            raise ValueError(f"Longitudinal family {family_id} has inconsistent n_ocms")
        primary_ids = [ocm_id for ocm_id in ids if ocm_id in by_ocm]
        missing_ids = [ocm_id for ocm_id in ids if ocm_id not in by_ocm]
        for ocm_id in primary_ids:
            previous = ocm_to_family.setdefault(ocm_id, family_id)
            if previous != family_id:
                raise ValueError(f"OCM {ocm_id} appears in multiple longitudinal families")
        primary_rows = [by_ocm[ocm_id] for ocm_id in primary_ids]
        outcome_ready = bool(primary_rows) and all(
            _boolean(row, "exact_paclitaxel_auc_available", context=row["run_accession"])
            and _boolean(row, "exact_paclitaxel_gi50_available", context=row["run_accession"])
            for row in primary_rows
        )
        exposure_ready = bool(primary_rows) and all(
            _boolean(
                row,
                "exact_cumulative_paclitaxel_exposure_available",
                context=row["run_accession"],
            )
            for row in primary_rows
        )
        if len(primary_ids) >= 2:
            pair_status = "pair_ready"
            analysis_status = "ready" if outcome_ready else "awaiting_exact_phenotype"
        elif len(primary_ids) == 1:
            pair_status = "partial_primary_family"
            analysis_status = "incomplete_primary_family"
        else:
            pair_status = "not_in_primary_cohort"
            analysis_status = "not_in_primary_cohort"
        design.append(
            {
                "family_id": family_id,
                "patient_id": family_id,
                "ocm_ids": ";".join(ids),
                "primary_ocm_ids": ";".join(primary_ids),
                "missing_primary_ocm_ids": ";".join(missing_ids),
                "primary_run_accessions": ";".join(
                    row["run_accession"] for row in primary_rows
                ),
                "n_family_ocms": str(len(ids)),
                "n_primary_ocms": str(len(primary_ids)),
                "relationship_type": family["relationship_type"],
                "analysis_role": family["analysis_role"],
                "pair_status": pair_status,
                "outcome_ready": str(outcome_ready).lower(),
                "exposure_ready": str(exposure_ready).lower(),
                "analysis_status": analysis_status,
            }
        )
    return design, ocm_to_family


def _track_status(
    rows: list[dict[str, str]],
    *,
    n_splits: int,
    outcome_fields: tuple[str, str],
) -> dict[str, Any]:
    patients = sorted({row["patient_id"] for row in rows})
    outcome_ready = bool(rows) and all(
        _boolean(row, field, context=row["run_accession"]) for row in rows for field in outcome_fields
    )
    result: dict[str, Any] = {
        "rows": len(rows),
        "ocms": len({row["canonical_ocm_id"] for row in rows}),
        "patients": len(patients),
        "outcome_ready": outcome_ready,
    }
    if not rows:
        result.update({"status": "blocked_missing_cohort_data", "cv": None})
        return result
    if len(patients) < n_splits:
        result.update({"status": "blocked_too_few_patient_groups", "cv": None})
        return result
    _assignment, loads = _assign_group_folds(rows, n_splits=n_splits)
    result["status"] = "ready" if outcome_ready else "awaiting_exact_phenotype"
    result["cv"] = {
        "method": "deterministic_greedy_group_assignment",
        "n_splits": n_splits,
        "fold_row_counts": loads,
        "group_column": "patient_id",
    }
    return result


def _write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str], *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temp_path = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp_path, path)


def _write_json(path: Path, value: dict[str, Any], *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        temp_path = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temp_path, path)


def prepare(
    *,
    manifest_path: Path,
    families_path: Path,
    schema_path: Path,
    output_dir: Path,
    overwrite: bool = False,
) -> dict[str, Any]:
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    rows = _read_tsv(manifest_path, REQUIRED_MANIFEST_FIELDS)
    expected_total = int(schema["cohort"]["all_runs"])
    if len(rows) != expected_total:
        raise ValueError(f"Expected {expected_total} manifest rows, found {len(rows)}")
    run_ids = [row["run_accession"] for row in rows]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("Manifest run_accession values are not unique")
    outcome_availability = _validate_outcomes(rows)
    primary = _primary_rows(rows, schema)
    families = _read_tsv(families_path, REQUIRED_FAMILY_FIELDS)
    family_design, ocm_to_family = _family_design(primary, families)

    n_splits = int(schema["validation"]["cv_splits"])
    cross_assignment, _cross_loads = _assign_group_folds(primary, n_splits=n_splits)
    intrinsic = []
    for row in primary:
        if not _boolean(
            row,
            "exact_cumulative_paclitaxel_exposure_available",
            context=row["run_accession"],
        ):
            continue
        exposure = _number(
            row,
            "cumulative_paclitaxel_exposure_mg",
            context=row["run_accession"],
            kind="nonnegative",
        )
        if exposure == 0:
            intrinsic.append(row)
    acquired_ocms = {
        ocm_id
        for family in family_design
        if family["pair_status"] == "pair_ready"
        for ocm_id in family["primary_ocm_ids"].split(";")
        if ocm_id
    }
    acquired = [row for row in primary if row["canonical_ocm_id"] in acquired_ocms]

    design_fields = [
        "study_accession",
        "run_accession",
        "canonical_ocm_id",
        "patient_id",
        "cross_sectional_eligible",
        "intrinsic_eligible",
        "acquired_family_id",
        "acquired_pair_eligible",
        "patient_cv_fold",
        "paclitaxel_auc_available",
        "paclitaxel_gi50_available",
        "cumulative_exposure_available",
    ]
    design_rows = []
    intrinsic_ids = {row["run_accession"] for row in intrinsic}
    for row in primary:
        design_rows.append(
            {
                "study_accession": row["study_accession"],
                "run_accession": row["run_accession"],
                "canonical_ocm_id": row["canonical_ocm_id"],
                "patient_id": row["patient_id"],
                "cross_sectional_eligible": "true",
                "intrinsic_eligible": str(row["run_accession"] in intrinsic_ids).lower(),
                "acquired_family_id": ocm_to_family.get(row["canonical_ocm_id"], ""),
                "acquired_pair_eligible": str(
                    row["canonical_ocm_id"] in acquired_ocms
                ).lower(),
                "patient_cv_fold": str(cross_assignment[row["patient_id"]]),
                "paclitaxel_auc_available": row["exact_paclitaxel_auc_available"],
                "paclitaxel_gi50_available": row["exact_paclitaxel_gi50_available"],
                "cumulative_exposure_available": row[
                    "exact_cumulative_paclitaxel_exposure_available"
                ],
            }
        )

    fold_fields = [
        "track",
        "patient_id",
        "fold",
        "n_primary_rows",
        "run_accessions",
    ]
    fold_rows: list[dict[str, str]] = []
    for track, track_rows in (
        ("cross_sectional", primary),
        ("intrinsic", intrinsic),
        ("acquired", acquired),
    ):
        if not track_rows:
            continue
        assignment, _loads = _assign_group_folds(track_rows, n_splits=n_splits)
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in track_rows:
            grouped[row["patient_id"]].append(row)
        for patient_id in sorted(grouped):
            members = grouped[patient_id]
            fold_rows.append(
                {
                    "track": track,
                    "patient_id": patient_id,
                    "fold": str(assignment[patient_id]),
                    "n_primary_rows": str(len(members)),
                    "run_accessions": ";".join(
                        sorted(row["run_accession"] for row in members)
                    ),
                }
            )

    outcome_fields = ("exact_paclitaxel_auc_available", "exact_paclitaxel_gi50_available")
    pair_ready = [family for family in family_design if family["pair_status"] == "pair_ready"]
    outcome_ready_pairs = [family for family in pair_ready if family["outcome_ready"] == "true"]
    summary = {
        "status": "prepared",
        "schema_version": schema["schema_version"],
        "manifest": {
            "path": str(manifest_path),
            "rows": len(rows),
            "study_counts": dict(sorted(Counter(row["study_accession"] for row in rows).items())),
            "unique_runs": len(set(run_ids)),
        },
        "primary_cohort": {
            "rows": len(primary),
            "unique_ocms": len({row["canonical_ocm_id"] for row in primary}),
            "unique_patients": len({row["patient_id"] for row in primary}),
            "sample_class": sorted({row["sample_class"] for row in primary}),
            "histotype_group": sorted({row["histotype_group"] for row in primary}),
        },
        "outcomes": outcome_availability,
        "tracks": {
            "cross_sectional": _track_status(
                primary, n_splits=n_splits, outcome_fields=outcome_fields
            ),
            "intrinsic": {
                **_track_status(intrinsic, n_splits=n_splits, outcome_fields=outcome_fields),
                "rule": schema["exposure"]["intrinsic_rule"],
                "proxy_not_used": schema["exposure"]["proxy_not_allowed"],
            },
            "acquired": {
                **_track_status(acquired, n_splits=n_splits, outcome_fields=outcome_fields),
                "longitudinal_families": len(family_design),
                "pair_ready_families": len(pair_ready),
                "outcome_ready_families": len(outcome_ready_pairs),
                "partial_or_absent_families": sum(
                    family["pair_status"] != "pair_ready" for family in family_design
                ),
            },
        },
        "cv": {
            "group_column": schema["validation"]["group_column"],
            "method": schema["validation"]["cv_method"],
            "n_splits": n_splits,
            "phenotype_used_for_assignment": False,
        },
        "outputs": {
            "design": str(output_dir / "taylor_design.tsv"),
            "patient_group_folds": str(output_dir / "patient_group_folds.tsv"),
            "acquired_families": str(output_dir / "acquired_family_design.tsv"),
            "summary": str(output_dir / "taylor_analysis_preparation.json"),
        },
    }

    _write_tsv(output_dir / "taylor_design.tsv", design_rows, design_fields, overwrite=overwrite)
    _write_tsv(
        output_dir / "patient_group_folds.tsv", fold_rows, fold_fields, overwrite=overwrite
    )
    family_fields = [
        "family_id",
        "patient_id",
        "ocm_ids",
        "primary_ocm_ids",
        "missing_primary_ocm_ids",
        "primary_run_accessions",
        "n_family_ocms",
        "n_primary_ocms",
        "relationship_type",
        "analysis_role",
        "pair_status",
        "outcome_ready",
        "exposure_ready",
        "analysis_status",
    ]
    _write_tsv(
        output_dir / "acquired_family_design.tsv",
        family_design,
        family_fields,
        overwrite=overwrite,
    )
    _write_json(
        output_dir / "taylor_analysis_preparation.json", summary, overwrite=overwrite
    )
    return summary


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "data/processed/metadata/ocm_master_manifest.tsv",
    )
    parser.add_argument(
        "--families",
        type=Path,
        default=root / "data/processed/metadata/longitudinal_families.tsv",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=root / "config/taylor_analysis.json",
    )
    parser.add_argument("--output-dir", type=Path, default=root / "data/processed/taylor")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    summary = prepare(
        manifest_path=args.manifest,
        families_path=args.families,
        schema_path=args.schema,
        output_dir=args.output_dir,
        overwrite=args.overwrite,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
