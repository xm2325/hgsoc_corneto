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
    product_paths = [Path(args.bgen), Path(args.sample), Path(args.bcftools_stats)]
    product_paths.extend(sorted(Path(args.parquet_dir).glob("*.parquet")))

    payload = {
        "source": {"url": args.source_url, "sha256": source_sha},
        "parameters": {"geno": args.geno, "maf": args.maf, "hwe": args.hwe},
        "tool_versions": {
            "bcftools": version(["bcftools", "--version"]),
            "plink2": version(["plink2", "--version"]),
            "python": version(["python", "--version"]),
        },
        "summary": summary,
        "products": {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in product_paths},
    }
    Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--source-url", required=True)
    p.add_argument("--source-sha-file", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--bgen", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--parquet-dir", required=True)
    p.add_argument("--bcftools-stats", required=True)
    p.add_argument("--geno", required=True)
    p.add_argument("--maf", required=True)
    p.add_argument("--hwe", required=True)
    p.add_argument("--output", required=True)
    return p


if __name__ == "__main__":
    main(parser().parse_args())
