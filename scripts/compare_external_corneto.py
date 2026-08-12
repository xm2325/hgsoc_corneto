#!/usr/bin/env python3
"""Fail-closed comparison of a frozen CORNETO signature across external groups.

This script deliberately consumes patient-level *network evidence*, not gene
expression.  Every external group must have been normalized and inferred
independently before it reaches this comparator.  A provenance contract and
complete patient-by-signature feature grid are mandatory, so raw TPM pooling
and treating cells as independent replicates cannot silently enter the result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


class ExternalValidationError(ValueError):
    """Raised when an external-validation input violates the frozen contract."""


ALLOWED_FEATURE_TYPES = {"edge", "reaction"}
ALLOWED_ANALYSIS_UNITS = {"patient_model", "patient_pseudobulk", "patient_tissue"}
ALLOWED_INPUT_SCALES = {
    "activity_score",
    "flux_direction",
    "network_selection",
    "standardized_rank",
}
SIGNATURE_COLUMNS = {"feature_type", "feature_id", "expected_direction"}
EVIDENCE_COLUMNS = {"patient_id", "feature_type", "feature_id", "selected", "direction"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path, required: set[str], label: str) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
    except OSError as error:
        raise ExternalValidationError(f"cannot read {label}: {error}") from error
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        missing = sorted(required - set(reader.fieldnames or []))
        raise ExternalValidationError(f"{label} is missing columns: {missing}")
    if not rows:
        raise ExternalValidationError(f"{label} is empty")
    return rows


def _direction(raw: str, label: str) -> int:
    try:
        value = int(raw)
    except (TypeError, ValueError) as error:
        raise ExternalValidationError(f"{label} must be -1, 0, or 1") from error
    if value not in {-1, 0, 1}:
        raise ExternalValidationError(f"{label} must be -1, 0, or 1")
    return value


def _boolean(raw: str, label: str) -> bool:
    value = str(raw).strip().casefold()
    if value in {"1", "true"}:
        return True
    if value in {"0", "false"}:
        return False
    raise ExternalValidationError(f"{label} must be 0/1 or true/false")


def read_signature(path: Path) -> dict[tuple[str, str], int]:
    rows = _read_tsv(path, SIGNATURE_COLUMNS, "frozen signature")
    result: dict[tuple[str, str], int] = {}
    for line, row in enumerate(rows, start=2):
        feature_type = row["feature_type"].strip().casefold()
        feature_id = row["feature_id"].strip()
        if feature_type not in ALLOWED_FEATURE_TYPES or not feature_id:
            raise ExternalValidationError(f"frozen signature line {line} has an invalid feature")
        key = (feature_type, feature_id)
        if key in result:
            raise ExternalValidationError(f"frozen signature has duplicate feature {key!r}")
        result[key] = _direction(
            row["expected_direction"], f"frozen signature line {line} expected_direction"
        )
    return result


def _object(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ExternalValidationError(f"{label} must be an object")
    return value


def _read_contract(
    path: Path,
    *,
    label: str,
    evidence_path: Path,
    signature_path: Path,
) -> dict[str, Any]:
    try:
        contract = _object(json.loads(path.read_text(encoding="utf-8")), f"{label} contract")
    except (OSError, json.JSONDecodeError) as error:
        raise ExternalValidationError(f"cannot read {label} contract: {error}") from error
    if contract.get("schema_version") != "external_corneto_group.v1":
        raise ExternalValidationError(f"{label}: unsupported contract schema")
    if contract.get("status") != "completed" or contract.get("group_label") != label:
        raise ExternalValidationError(f"{label}: contract is not completed or is mislabeled")
    if not isinstance(contract.get("source_accession"), str) or not contract["source_accession"]:
        raise ExternalValidationError(f"{label}: source_accession is missing")
    if contract.get("evidence_sha256") != _sha256(evidence_path):
        raise ExternalValidationError(f"{label}: evidence SHA-256 mismatch")
    if contract.get("signature_sha256") != _sha256(signature_path):
        raise ExternalValidationError(f"{label}: frozen signature SHA-256 mismatch")
    if contract.get("analysis_unit") not in ALLOWED_ANALYSIS_UNITS:
        raise ExternalValidationError(f"{label}: analysis_unit must be patient-level")

    normalization = _object(contract.get("normalization"), f"{label}.normalization")
    if normalization.get("performed_within_dataset") is not True:
        raise ExternalValidationError(f"{label}: normalization was not performed within dataset")
    if normalization.get("pooled_raw_expression") is not False:
        raise ExternalValidationError(f"{label}: pooled raw expression is forbidden")
    if normalization.get("input_scale") not in ALLOWED_INPUT_SCALES:
        raise ExternalValidationError(f"{label}: unsupported normalized input scale")

    independence = _object(contract.get("independence"), f"{label}.independence")
    if independence.get("cells_as_replicates") is not False:
        raise ExternalValidationError(f"{label}: cells-as-n is forbidden")
    if independence.get("patient_id_column") != "patient_id":
        raise ExternalValidationError(f"{label}: patient_id_column must be 'patient_id'")

    inference = _object(contract.get("inference"), f"{label}.inference")
    if inference.get("signature_frozen_before_external_scoring") is not True:
        raise ExternalValidationError(f"{label}: signature was not frozen before scoring")
    if inference.get("feature_selection_using_external_labels") is not False:
        raise ExternalValidationError(f"{label}: external-label feature selection is forbidden")
    return contract


def _read_group(
    label: str,
    evidence_path: Path,
    contract_path: Path,
    signature_path: Path,
    signature: dict[tuple[str, str], int],
    min_patients: int,
) -> dict[str, Any]:
    contract = _read_contract(
        contract_path,
        label=label,
        evidence_path=evidence_path,
        signature_path=signature_path,
    )
    rows = _read_tsv(evidence_path, EVIDENCE_COLUMNS, f"{label} evidence")
    patient_rows: dict[str, dict[tuple[str, str], tuple[bool, int]]] = {}
    seen: set[tuple[str, str, str]] = set()
    for line, row in enumerate(rows, start=2):
        patient = row["patient_id"].strip()
        feature_type = row["feature_type"].strip().casefold()
        feature_id = row["feature_id"].strip()
        feature = (feature_type, feature_id)
        if not patient:
            raise ExternalValidationError(f"{label} evidence line {line} has no patient_id")
        if feature not in signature:
            raise ExternalValidationError(
                f"{label} evidence line {line} contains non-frozen feature {feature!r}"
            )
        key = (patient, feature_type, feature_id)
        if key in seen:
            raise ExternalValidationError(f"{label} evidence duplicates patient-feature {key!r}")
        seen.add(key)
        selected = _boolean(row["selected"], f"{label} evidence line {line} selected")
        direction = _direction(row["direction"], f"{label} evidence line {line} direction")
        if not selected and direction != 0:
            raise ExternalValidationError(
                f"{label} evidence line {line}: unselected feature must have direction 0"
            )
        if selected and direction == 0 and signature[feature] != 0:
            raise ExternalValidationError(
                f"{label} evidence line {line}: selected directed feature has direction 0"
            )
        patient_rows.setdefault(patient, {})[feature] = (selected, direction)

    patients = sorted(patient_rows)
    if len(patients) < min_patients:
        raise ExternalValidationError(
            f"{label}: {len(patients)} patients, fewer than required {min_patients}"
        )
    expected = set(signature)
    for patient in patients:
        observed = set(patient_rows[patient])
        if observed != expected:
            missing = sorted(expected - observed)[:5]
            raise ExternalValidationError(
                f"{label}: patient {patient!r} does not have the complete frozen grid; "
                f"missing examples={missing}, extra_count={len(observed - expected)}"
            )
    if contract.get("patient_count") != len(patients):
        raise ExternalValidationError(f"{label}: contract patient_count disagrees with evidence")
    return {
        "label": label,
        "evidence_path": evidence_path,
        "contract_path": contract_path,
        "contract": contract,
        "patients": patients,
        "rows": patient_rows,
    }


def _quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ExternalValidationError("cannot take quantile of an empty sequence")
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def _jaccard(left: set[tuple[str, str]], right: set[tuple[str, str]]) -> float | None:
    union = left | right
    return None if not union else len(left & right) / len(union)


def _summarize_draw(
    rows: dict[str, dict[tuple[str, str], tuple[bool, int]]],
    sampled_patients: list[str],
    features: list[tuple[str, str]],
) -> dict[tuple[str, str], tuple[float, int, int]]:
    result: dict[tuple[str, str], tuple[float, int, int]] = {}
    denominator = len(sampled_patients)
    for feature in features:
        values = [rows[patient][feature] for patient in sampled_patients]
        selected_count = sum(selected for selected, _ in values)
        signed_sum = sum(direction for selected, direction in values if selected)
        consensus_direction = 1 if signed_sum > 0 else -1 if signed_sum < 0 else 0
        result[feature] = (selected_count / denominator, consensus_direction, selected_count)
    return result


def compare(
    *,
    signature_path: Path,
    group_specs: list[str],
    bootstrap_iterations: int = 1000,
    seed: int = 1729,
    prevalence_threshold: float = 0.5,
    min_patients: int = 2,
) -> dict[str, Any]:
    if bootstrap_iterations < 100:
        raise ExternalValidationError("bootstrap_iterations must be at least 100")
    if not 0 < prevalence_threshold <= 1:
        raise ExternalValidationError("prevalence_threshold must be in (0, 1]")
    if min_patients < 2:
        raise ExternalValidationError("min_patients must be at least 2")
    signature = read_signature(signature_path)
    parsed: list[tuple[str, Path, Path]] = []
    for spec in group_specs:
        parts = spec.split("=", 1)
        if len(parts) != 2 or "," not in parts[1]:
            raise ExternalValidationError(
                f"group must be LABEL=EVIDENCE_TSV,CONTRACT_JSON, got {spec!r}"
            )
        label = parts[0].strip()
        evidence, contract = parts[1].split(",", 1)
        if not label or not evidence or not contract:
            raise ExternalValidationError(f"invalid group specification {spec!r}")
        parsed.append((label, Path(evidence), Path(contract)))
    labels = [item[0] for item in parsed]
    if len(labels) < 2 or len(labels) != len(set(labels)):
        raise ExternalValidationError("at least two uniquely labeled groups are required")
    groups = {
        label: _read_group(
            label,
            evidence,
            contract,
            signature_path,
            signature,
            min_patients,
        )
        for label, evidence, contract in parsed
    }
    features = sorted(signature)
    rng = random.Random(seed)
    bootstrap: dict[str, list[dict[tuple[str, str], tuple[float, int, int]]]] = {}
    summaries: dict[str, Any] = {}
    for label, group in groups.items():
        observed = _summarize_draw(group["rows"], group["patients"], features)
        draws = []
        for _ in range(bootstrap_iterations):
            sampled = [rng.choice(group["patients"]) for _ in group["patients"]]
            draws.append(_summarize_draw(group["rows"], sampled, features))
        bootstrap[label] = draws
        feature_rows = []
        consensus: set[tuple[str, str]] = set()
        for feature in features:
            prevalence, consensus_direction, selected_count = observed[feature]
            if prevalence >= prevalence_threshold:
                consensus.add(feature)
            direction_observations = [
                direction
                for patient in group["patients"]
                for selected, direction in [group["rows"][patient][feature]]
                if selected and direction != 0
            ]
            expected_direction = signature[feature]
            direction_concordance = None
            if expected_direction != 0 and direction_observations:
                direction_concordance = sum(
                    value == expected_direction for value in direction_observations
                ) / len(direction_observations)
            prevalence_draws = [draw[feature][0] for draw in draws]
            feature_rows.append(
                {
                    "feature_type": feature[0],
                    "feature_id": feature[1],
                    "expected_direction": expected_direction,
                    "patient_count": len(group["patients"]),
                    "selected_patient_count": selected_count,
                    "prevalence": prevalence,
                    "prevalence_ci95": [
                        _quantile(prevalence_draws, 0.025),
                        _quantile(prevalence_draws, 0.975),
                    ],
                    "consensus_direction": consensus_direction,
                    "direction_concordance": direction_concordance,
                }
            )
        summaries[label] = {
            "source_accession": group["contract"]["source_accession"],
            "analysis_unit": group["contract"]["analysis_unit"],
            "patient_count": len(group["patients"]),
            "consensus_feature_count": len(consensus),
            "signature_jaccard": _jaccard(consensus, set(features)),
            "features": feature_rows,
            "_consensus": consensus,
            "_observed": observed,
        }

    pairwise = []
    for index, left in enumerate(labels):
        for right in labels[index + 1 :]:
            left_consensus = summaries[left]["_consensus"]
            right_consensus = summaries[right]["_consensus"]
            shared = left_consensus & right_consensus
            directed_shared = [
                feature
                for feature in shared
                if summaries[left]["_observed"][feature][1] != 0
                and summaries[right]["_observed"][feature][1] != 0
            ]
            direction_concordance = None
            if directed_shared:
                direction_concordance = sum(
                    summaries[left]["_observed"][feature][1]
                    == summaries[right]["_observed"][feature][1]
                    for feature in directed_shared
                ) / len(directed_shared)
            jaccard_draws = []
            for left_draw, right_draw in zip(bootstrap[left], bootstrap[right], strict=True):
                left_set = {
                    feature for feature in features if left_draw[feature][0] >= prevalence_threshold
                }
                right_set = {
                    feature for feature in features if right_draw[feature][0] >= prevalence_threshold
                }
                value = _jaccard(left_set, right_set)
                if value is not None:
                    jaccard_draws.append(value)

            common_patients = sorted(set(groups[left]["patients"]) & set(groups[right]["patients"]))
            patient_jaccards = []
            for patient in common_patients:
                left_set = {
                    feature for feature in features if groups[left]["rows"][patient][feature][0]
                }
                right_set = {
                    feature for feature in features if groups[right]["rows"][patient][feature][0]
                }
                value = _jaccard(left_set, right_set)
                if value is not None:
                    patient_jaccards.append(value)
            pairwise.append(
                {
                    "left": left,
                    "right": right,
                    "consensus_intersection": len(shared),
                    "consensus_union": len(left_consensus | right_consensus),
                    "consensus_jaccard": _jaccard(left_consensus, right_consensus),
                    "consensus_jaccard_ci95": (
                        [_quantile(jaccard_draws, 0.025), _quantile(jaccard_draws, 0.975)]
                        if jaccard_draws
                        else None
                    ),
                    "shared_directed_feature_count": len(directed_shared),
                    "direction_concordance": direction_concordance,
                    "matched_patient_count": len(common_patients),
                    "mean_matched_patient_jaccard": (
                        sum(patient_jaccards) / len(patient_jaccards)
                        if patient_jaccards
                        else None
                    ),
                }
            )

    for summary in summaries.values():
        del summary["_consensus"]
        del summary["_observed"]
    return {
        "schema_version": "external_corneto_comparison.v1",
        "status": "completed",
        "completed_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
        "frozen_signature": {
            "path": str(signature_path),
            "sha256": _sha256(signature_path),
            "feature_count": len(signature),
            "feature_type_counts": {
                feature_type: sum(feature[0] == feature_type for feature in signature)
                for feature_type in sorted({feature[0] for feature in signature})
            },
        },
        "contract": {
            "normalization": "independent within each source dataset before CORNETO scoring",
            "analysis_unit": "patient; single cells may only contribute through patient pseudobulk",
            "raw_expression_pooling": "forbidden",
            "signature_selection": "frozen before external scoring",
            "prevalence_threshold": prevalence_threshold,
            "bootstrap_unit": "patient_id",
            "bootstrap_iterations": bootstrap_iterations,
            "random_seed": seed,
        },
        "groups": summaries,
        "pairwise": pairwise,
        "claim_limit": (
            "External concordance of a prespecified CORNETO signature. This is not measured flux, "
            "causal evidence, drug-response validation, or permission to treat cells as replicates."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument(
        "--group",
        action="append",
        required=True,
        metavar="LABEL=EVIDENCE_TSV,CONTRACT_JSON",
    )
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--prevalence-threshold", type=float, default=0.5)
    parser.add_argument("--min-patients", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        print(f"refusing to overwrite output: {args.output}", file=sys.stderr)
        return 1
    try:
        result = compare(
            signature_path=args.signature,
            group_specs=args.group,
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
            prevalence_threshold=args.prevalence_threshold,
            min_patients=args.min_patients,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    except (ExternalValidationError, OSError) as error:
        print(json.dumps({"status": "invalid", "error": str(error)}), file=sys.stderr)
        return 1
    print(json.dumps({"status": "completed", "output": str(args.output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
