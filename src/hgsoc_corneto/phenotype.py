"""Fail-closed validation for an externally supplied Taylor phenotype table.

This module deliberately does not fetch, infer, impute, transform, merge, or
analyse phenotype values.  It validates an exact, source-attributed intake
table against the frozen primary-cohort OCM-to-patient mapping and returns a
value-free readiness receipt.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

MISSING = frozenset({"", "NA", "N/A", "null", "None"})
REQUIRED_INTAKE_FIELDS = frozenset(
    {
        "canonical_ocm_id",
        "patient_id",
        "drug",
        "endpoint",
        "value",
        "unit",
        "endpoint_definition",
        "is_exact",
        "source_file",
        "source_record_id",
    }
)
REQUIRED_MANIFEST_FIELDS = frozenset(
    {
        "run_accession",
        "canonical_ocm_id",
        "patient_id",
        "primary_cohort_eligible",
    }
)


class PhenotypeIntakeError(ValueError):
    """Raised when the supplied phenotype table is unsafe to admit."""


def _read_tsv(path: Path, required: frozenset[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PhenotypeIntakeError(f"Missing TSV header: {path}")
        fields = list(reader.fieldnames)
        if len(fields) != len(set(fields)):
            raise PhenotypeIntakeError(f"Duplicate TSV columns: {path}")
        missing = sorted(required - set(fields))
        if missing:
            raise PhenotypeIntakeError(f"Missing required columns in {path}: {missing}")
        return list(reader)


def _bool(value: str, *, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise PhenotypeIntakeError(f"{context}: expected true/false, got {value!r}")


def _number(value: str, *, context: str, constraint: str) -> float:
    if value.strip() in MISSING:
        raise PhenotypeIntakeError(f"{context}: exact value is missing")
    try:
        parsed = float(value)
    except ValueError as error:
        raise PhenotypeIntakeError(f"{context}: value is not numeric: {value!r}") from error
    if not math.isfinite(parsed):
        raise PhenotypeIntakeError(f"{context}: value must be finite")
    if constraint == "nonnegative" and parsed < 0:
        raise PhenotypeIntakeError(f"{context}: value must be non-negative")
    if constraint == "positive" and parsed <= 0:
        raise PhenotypeIntakeError(f"{context}: value must be positive")
    if constraint not in {"nonnegative", "positive"}:
        raise PhenotypeIntakeError(f"Unsupported numeric constraint: {constraint!r}")
    return parsed


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _primary_mapping(manifest_path: Path, expected_ocms: int) -> dict[str, str]:
    rows = _read_tsv(manifest_path, REQUIRED_MANIFEST_FIELDS)
    primary: list[dict[str, str]] = []
    for row in rows:
        context = row.get("run_accession", "manifest row")
        if _bool(row["primary_cohort_eligible"], context=context):
            primary.append(row)
    if len(primary) != expected_ocms:
        raise PhenotypeIntakeError(
            f"Frozen cohort expects {expected_ocms} primary rows, found {len(primary)}"
        )
    mapping: dict[str, str] = {}
    for row in primary:
        ocm_id = row["canonical_ocm_id"].strip()
        patient_id = row["patient_id"].strip()
        if not ocm_id or not patient_id:
            raise PhenotypeIntakeError("Primary manifest row has missing OCM/patient mapping")
        if ocm_id in mapping:
            raise PhenotypeIntakeError(f"Duplicate primary OCM in manifest: {ocm_id}")
        mapping[ocm_id] = patient_id
    return mapping


def blocked_receipt(
    *, phenotype_path: Path, manifest_path: Path, schema_path: Path
) -> dict[str, Any]:
    """Return a value-free receipt when the external table does not exist."""

    return {
        "schema_version": 1,
        "status": "blocked_missing_phenotype_file",
        "association_allowed": False,
        "reason": "Exact Taylor phenotype intake was not supplied; no values were inferred.",
        "inputs": {
            "phenotype": str(phenotype_path),
            "manifest": str(manifest_path),
            "schema": str(schema_path),
        },
        "phenotype_values_written": False,
        "association_run": False,
    }


def validate_phenotype_intake(
    *, phenotype_path: Path, manifest_path: Path, schema_path: Path
) -> dict[str, Any]:
    """Validate exact phenotype rows and return a receipt containing no values.

    Missing endpoint rows produce a valid but blocked receipt.  Unsafe content
    (duplicates, unit mismatches, non-exact rows, unknown OCMs, or patient-map
    conflicts) raises :class:`PhenotypeIntakeError`.
    """

    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    endpoints: dict[str, dict[str, Any]] = schema["endpoints"]
    expected_ocms = int(schema["cohort"]["primary_ocms"])
    mapping = _primary_mapping(manifest_path, expected_ocms)
    rows = _read_tsv(phenotype_path, REQUIRED_INTAKE_FIELDS)

    seen: set[tuple[str, str]] = set()
    seen_source_records: set[tuple[str, str, str, str]] = set()
    endpoint_ocms: dict[str, set[str]] = defaultdict(set)
    definitions: dict[str, set[str]] = defaultdict(set)

    for index, row in enumerate(rows, start=2):
        endpoint = row["endpoint"].strip()
        context = f"{phenotype_path.name}:{index}"
        if endpoint not in endpoints:
            raise PhenotypeIntakeError(f"{context}: unsupported endpoint {endpoint!r}")
        spec = endpoints[endpoint]
        if row["drug"].strip().lower() != str(spec["drug"]).lower():
            raise PhenotypeIntakeError(f"{context}: endpoint is not labelled as paclitaxel")
        if not _bool(row["is_exact"], context=context):
            raise PhenotypeIntakeError(f"{context}: proxy/imputed rows are forbidden")

        ocm_id = row["canonical_ocm_id"].strip()
        patient_id = row["patient_id"].strip()
        if ocm_id not in mapping:
            raise PhenotypeIntakeError(
                f"{context}: OCM {ocm_id!r} is not in the frozen primary cohort"
            )
        if patient_id != mapping[ocm_id]:
            raise PhenotypeIntakeError(
                f"{context}: patient mapping conflict for {ocm_id}: "
                f"expected {mapping[ocm_id]!r}, got {patient_id!r}"
            )

        key = (ocm_id, endpoint)
        if key in seen:
            raise PhenotypeIntakeError(f"{context}: duplicate OCM/endpoint row: {key}")
        seen.add(key)
        source_key = (
            row["source_file"].strip(),
            row["source_record_id"].strip(),
            ocm_id,
            endpoint,
        )
        if any(value in MISSING for value in source_key[:2]):
            raise PhenotypeIntakeError(f"{context}: source provenance is missing")
        if source_key in seen_source_records:
            raise PhenotypeIntakeError(f"{context}: duplicate source record")
        seen_source_records.add(source_key)

        unit = row["unit"].strip()
        accepted_units = {str(value) for value in spec["accepted_units"]}
        if not accepted_units:
            raise PhenotypeIntakeError(
                f"{endpoint}: schema has no source-confirmed accepted unit; update schema first"
            )
        if unit not in accepted_units:
            raise PhenotypeIntakeError(
                f"{context}: unit {unit!r} is not accepted for {endpoint}; "
                f"expected one of {sorted(accepted_units)}"
            )
        definition = row["endpoint_definition"].strip()
        if definition in MISSING:
            raise PhenotypeIntakeError(f"{context}: endpoint definition is missing")

        _number(row["value"], context=context, constraint=str(spec["constraint"]))
        endpoint_ocms[endpoint].add(ocm_id)
        definitions[endpoint].add(definition)

    endpoint_summary: dict[str, dict[str, Any]] = {}
    for endpoint, spec in endpoints.items():
        observed = endpoint_ocms.get(endpoint, set())
        missing = sorted(set(mapping) - observed)
        endpoint_summary[endpoint] = {
            "role": spec["role"],
            "transform_after_gate": spec.get("transform_after_gate"),
            "rows": len(observed),
            "expected_rows": expected_ocms,
            "complete": not missing,
            "missing_ocm_ids": missing,
            "units": sorted(set(spec["accepted_units"])),
            "endpoint_definition_count": len(definitions.get(endpoint, set())),
        }

    primary_ready = endpoint_summary[str(schema["analysis_readiness"]["primary_endpoint"])][
        "complete"
    ]
    secondary_ready = endpoint_summary[
        str(schema["analysis_readiness"]["secondary_endpoint"])
    ]["complete"]
    exposure_ready = endpoint_summary[str(schema["analysis_readiness"]["exposure_endpoint"])][
        "complete"
    ]
    all_ready = primary_ready and secondary_ready and exposure_ready
    return {
        "schema_version": schema["schema_version"],
        "status": "ready" if all_ready else "blocked_incomplete_exact_phenotype",
        "association_allowed": bool(primary_ready),
        "secondary_analysis_allowed": bool(primary_ready and secondary_ready),
        "exposure_stratified_analysis_allowed": bool(primary_ready and exposure_ready),
        "all_endpoints_ready": bool(all_ready),
        "cohort": {
            "expected_primary_ocms": expected_ocms,
            "mapped_primary_ocms": len(mapping),
        },
        "endpoints": endpoint_summary,
        "input_fingerprint": {
            "path": str(phenotype_path),
            "sha256": _sha256(phenotype_path),
            "rows": len(rows),
            "endpoint_counts": dict(sorted(Counter(row["endpoint"] for row in rows).items())),
        },
        "safety": {
            "patient_mapping_validated": True,
            "duplicates_rejected": True,
            "units_validated": True,
            "exact_only": True,
            "phenotype_values_written": False,
            "association_run": False,
        },
    }
