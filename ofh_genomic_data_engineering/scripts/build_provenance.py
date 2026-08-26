from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        text = (result.stdout or result.stderr).strip().splitlines()
        return text[0] if text else "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def main(args: argparse.Namespace) -> None:
    source_sha = Path(args.source_sha_file).read_text().split()[0]
    summary = json.loads(Path(args.summary).read_text())
    delivery_validation = json.loads(Path(args.delivery_validation).read_text())
    bgen_validation = json.loads(Path(args.bgen_validation).read_text())
    metadata_validation = json.loads(Path(args.metadata_validation).read_text())
    delivery = delivery_validation.get("delivery", {})
    observed = delivery_validation.get("source_observed", {})

    if delivery_validation.get("status") != "PASS":
        raise ValueError("provenance cannot be built from a failed delivery validation")
    if delivery_validation.get("action") != "PROCESS":
        raise ValueError("provenance requires a PROCESS delivery decision")
    if observed.get("sha256") != source_sha:
        raise ValueError("delivery source SHA does not match downloaded source SHA")
    if bgen_validation.get("status") != "PASS":
        raise ValueError("provenance cannot be built from a failed BGEN round-trip validation")
    if metadata_validation.get("status") != "PASS":
        raise ValueError("provenance cannot be built from a failed sample metadata validation")

    product_paths = [
        Path(args.bgen),
        Path(args.sample),
        Path(args.bgen_validation),
        Path(args.sample_metadata),
        Path(args.metadata_validation),
        Path(args.bcftools_stats),
        Path(args.schema_manifest),
        Path(args.query_validation),
        Path(args.delivery_validation),
    ]
    product_paths.extend(sorted(Path(args.parquet_dir).glob("*.parquet")))

    payload = {
        "source": {"url": args.source_url, "sha256": source_sha},
        "delivery": {
            "delivery_id": delivery.get("delivery_id"),
            "provider": delivery.get("provider"),
            "reference_genome": delivery.get("reference_genome"),
            "delivery_fingerprint": delivery.get("delivery_fingerprint"),
            "sample_count": observed.get("sample_count"),
            "sample_ids_sha256": observed.get("sample_ids_sha256"),
        },
        "parameters": {
            "geno": args.geno,
            "maf": args.maf,
            "hwe": args.hwe,
            "plink_seed": int(args.plink_seed),
            "plink_threads": int(args.plink_threads),
            "plink_memory_mb": int(args.plink_memory_mb),
            "delivery_fingerprint": delivery.get("delivery_fingerprint"),
            "reference_genome": delivery.get("reference_genome"),
        },
        "bgen_roundtrip": bgen_validation,
        "sample_metadata_join": metadata_validation,
        "tool_versions": {
            "bcftools": version(["bcftools", "--version"]),
            "plink2": version(["plink2", "--version"]),
            "python": version(["python", "--version"]),
            "duckdb": version(["python", "-c", "import duckdb; print(duckdb.__version__)"]),
        },
        "summary": summary,
        "products": {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in product_paths},
    }
    Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--source-url", required=True)
    p.add_argument("--source-sha-file", required=True)
    p.add_argument("--delivery-validation", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--bgen", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--bgen-validation", required=True)
    p.add_argument("--sample-metadata", required=True)
    p.add_argument("--metadata-validation", required=True)
    p.add_argument("--parquet-dir", required=True)
    p.add_argument("--bcftools-stats", required=True)
    p.add_argument("--schema-manifest", required=True)
    p.add_argument("--query-validation", required=True)
    p.add_argument("--geno", required=True)
    p.add_argument("--maf", required=True)
    p.add_argument("--hwe", required=True)
    p.add_argument("--plink-seed", required=True)
    p.add_argument("--plink-threads", required=True)
    p.add_argument("--plink-memory-mb", required=True)
    p.add_argument("--output", required=True)
    return p


if __name__ == "__main__":
    main(parser().parse_args())
