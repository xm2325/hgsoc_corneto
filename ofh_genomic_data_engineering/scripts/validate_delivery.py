from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any


REQUIRED_FIELDS = (
    "schema_version",
    "delivery_id",
    "provider",
    "source_uri",
    "source_format",
    "reference_genome",
    "sha256",
    "sample_roster",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_sha256(payload: object) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def sample_ids_sha256(sample_ids: list[str]) -> str:
    raw = ("\n".join(sorted(sample_ids)) + "\n").encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def read_vcf_sample_ids(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                if len(columns) < 9:
                    raise ValueError("VCF header has fewer than 9 fixed columns")
                return columns[9:]
    raise ValueError("VCF #CHROM header not found")


def load_registry(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    payload = json.loads(path.read_text())
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict) and isinstance(payload.get("deliveries"), list):
        return payload["deliveries"]
    raise ValueError("registry must be a list or an object with a deliveries list")


def validate_delivery(
    *,
    manifest_path: Path,
    source_path: Path,
    output_path: Path,
    required_reference_genome: str,
    expected_source_uri: str | None = None,
    registry_path: Path | None = None,
) -> dict[str, object]:
    manifest = json.loads(manifest_path.read_text())
    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, detail: object) -> None:
        checks.append({"name": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    missing_fields = [field for field in REQUIRED_FIELDS if field not in manifest]
    record("manifest.required_fields", not missing_fields, {"missing": missing_fields})

    schema_version = manifest.get("schema_version")
    record("manifest.schema_version", schema_version == 1, schema_version)

    delivery_id = manifest.get("delivery_id")
    record("manifest.delivery_id", isinstance(delivery_id, str) and bool(delivery_id.strip()), delivery_id)

    provider = manifest.get("provider")
    record("manifest.provider", isinstance(provider, str) and bool(provider.strip()), provider)

    source_uri = manifest.get("source_uri")
    record("manifest.source_uri", isinstance(source_uri, str) and bool(source_uri.strip()), source_uri)
    if expected_source_uri is not None:
        record(
            "manifest.source_uri_matches_pipeline",
            source_uri == expected_source_uri,
            {"manifest": source_uri, "pipeline": expected_source_uri},
        )

    source_format = manifest.get("source_format")
    record("manifest.source_format", source_format == "vcf.gz", source_format)

    reference_genome = manifest.get("reference_genome")
    record(
        "manifest.reference_genome",
        reference_genome == required_reference_genome,
        {"declared": reference_genome, "required": required_reference_genome},
    )

    source_exists = source_path.is_file() and source_path.stat().st_size > 0
    record("source.file_nonempty", source_exists, str(source_path))

    observed_sha = sha256(source_path) if source_exists else None
    expected_sha = manifest.get("sha256")
    record(
        "source.sha256",
        isinstance(expected_sha, str) and len(expected_sha) == 64 and observed_sha == expected_sha,
        {"expected": expected_sha, "observed": observed_sha},
    )

    sample_ids: list[str] = []
    sample_error: str | None = None
    if source_exists:
        try:
            sample_ids = read_vcf_sample_ids(source_path)
        except (OSError, EOFError, ValueError) as exc:
            sample_error = str(exc)
    record(
        "source.vcf_header",
        sample_error is None and bool(sample_ids),
        {"error": sample_error, "sample_count": len(sample_ids)},
    )

    duplicate_samples = len(sample_ids) - len(set(sample_ids))
    record(
        "source.sample_ids_unique",
        bool(sample_ids) and duplicate_samples == 0,
        {"sample_count": len(sample_ids), "duplicate_count": duplicate_samples},
    )

    roster = manifest.get("sample_roster")
    roster = roster if isinstance(roster, dict) else {}
    expected_count = roster.get("count")
    observed_roster_sha = sample_ids_sha256(sample_ids) if sample_ids else None
    expected_roster_sha = roster.get("ids_sha256")
    record(
        "source.sample_count",
        isinstance(expected_count, int) and expected_count == len(sample_ids),
        {"expected": expected_count, "observed": len(sample_ids)},
    )
    record(
        "source.sample_roster_sha256",
        isinstance(expected_roster_sha, str)
        and len(expected_roster_sha) == 64
        and observed_roster_sha == expected_roster_sha,
        {"expected": expected_roster_sha, "observed": observed_roster_sha},
    )

    delivery_basis = {
        "schema_version": schema_version,
        "delivery_id": delivery_id,
        "provider": provider,
        "source_uri": source_uri,
        "source_format": source_format,
        "reference_genome": reference_genome,
        "source_sha256": observed_sha,
        "sample_roster": {"count": len(sample_ids), "ids_sha256": observed_roster_sha},
    }
    delivery_fingerprint = canonical_sha256(delivery_basis)

    registry_entries: list[dict[str, Any]] = []
    registry_error: str | None = None
    try:
        registry_entries = load_registry(registry_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        registry_error = str(exc)
    record("registry.readable", registry_error is None, {"error": registry_error})

    matching = [
        item for item in registry_entries
        if isinstance(item, dict) and item.get("delivery_id") == delivery_id
    ]
    record(
        "registry.delivery_id_unique",
        len(matching) <= 1,
        {"delivery_id": delivery_id, "matches": len(matching)},
    )

    action = "PROCESS"
    should_process = True
    if len(matching) == 1:
        prior = matching[0]
        same_content = (
            prior.get("source_sha256") == observed_sha
            and prior.get("delivery_fingerprint") == delivery_fingerprint
        )
        record(
            "registry.delivery_id_collision",
            same_content,
            {
                "delivery_id": delivery_id,
                "same_content": same_content,
                "registered_source_sha256": prior.get("source_sha256"),
                "observed_source_sha256": observed_sha,
            },
        )
        if same_content:
            action = "NOOP"
            should_process = False
        else:
            action = "REJECT"
            should_process = False
    elif len(matching) > 1:
        action = "REJECT"
        should_process = False
    else:
        record(
            "registry.delivery_id_collision",
            True,
            {"delivery_id": delivery_id, "same_content": None},
        )

    passed = all(check["status"] == "PASS" for check in checks)
    if not passed:
        action = "REJECT"
        should_process = False

    payload = {
        "status": "PASS" if passed else "FAIL",
        "action": action,
        "should_process": should_process,
        "delivery": {
            "delivery_id": delivery_id,
            "provider": provider,
            "source_uri": source_uri,
            "source_format": source_format,
            "reference_genome": reference_genome,
            "delivery_fingerprint": delivery_fingerprint,
        },
        "source_observed": {
            "sha256": observed_sha,
            "sample_count": len(sample_ids),
            "sample_ids_sha256": observed_roster_sha,
        },
        "checks": checks,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--source", required=True)
    p.add_argument("--required-reference-genome", required=True)
    p.add_argument("--expected-source-uri")
    p.add_argument("--registry")
    p.add_argument("--output", required=True)
    return p


def main() -> int:
    args = parser().parse_args()
    payload = validate_delivery(
        manifest_path=Path(args.manifest),
        source_path=Path(args.source),
        output_path=Path(args.output),
        required_reference_genome=args.required_reference_genome,
        expected_source_uri=args.expected_source_uri,
        registry_path=Path(args.registry) if args.registry else None,
    )
    if payload["status"] != "PASS":
        failed = [check["name"] for check in payload["checks"] if check["status"] == "FAIL"]
        print("delivery validation failed: " + ", ".join(failed))
        return 2
    print(
        "delivery validation passed: "
        f"{payload['action']} {payload['delivery']['delivery_fingerprint']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
