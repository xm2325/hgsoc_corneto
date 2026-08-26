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


def run(command: list[str]) -> str:
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    return result.stdout.strip()


def inspect_vcf(vcf: Path, bcftools: str = "bcftools") -> dict[str, object]:
    sample_ids = [line for line in run([bcftools, "query", "-l", str(vcf)]).splitlines() if line]
    if not sample_ids:
        raise ValueError("input VCF contains no samples")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("input VCF contains duplicate sample IDs")

    variant_count_text = run([bcftools, "index", "-n", str(vcf)])
    variant_count = int(variant_count_text)
    if variant_count <= 0:
        raise ValueError("input VCF contains no indexed variants")

    sample_payload = ("\n".join(sample_ids) + "\n").encode("utf-8")
    return {
        "sample_count": len(sample_ids),
        "sample_ids_sha256": hashlib.sha256(sample_payload).hexdigest(),
        "variant_count": variant_count,
        "normalised_vcf_sha256": sha256(vcf),
    }


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--bcftools", default="bcftools")
    return p


def main() -> None:
    args = parser().parse_args()
    payload = inspect_vcf(Path(args.vcf), args.bcftools)
    Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
