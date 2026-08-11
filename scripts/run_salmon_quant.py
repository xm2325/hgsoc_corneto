#!/usr/bin/env python3
"""Download, verify, and quantify one ENA paired-end RNA-seq run."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from hgsoc_corneto.io import write_json
from hgsoc_corneto.rna import (
    FastqSpec,
    load_rna_run_specs,
    resumable_curl_command,
    validate_fastq_file,
)


def _run(command: list[str], *, capture: bool = False) -> str:
    result = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
    )
    return result.stdout.strip() if capture else ""


def _download_fastq(spec: FastqSpec, target: Path) -> dict[str, Any]:
    valid, reason = validate_fastq_file(target, spec)
    if not valid and target.exists():
        raise ValueError(f"Existing FASTQ requires inspection: {target} ({reason})")
    if not valid:
        partial = target.with_name(target.name + ".partial")
        _run(resumable_curl_command(url=spec.url, target=partial))
        valid, reason = validate_fastq_file(partial, spec)
        if not valid:
            raise ValueError(f"Downloaded FASTQ failed verification: {partial} ({reason})")
        os.replace(partial, target)
    return {
        "mate": spec.mate,
        "url": spec.url,
        "path": str(target),
        "bytes": spec.bytes,
        "md5": spec.md5,
        "verification": "verified",
    }


def _validate_existing_quant(output: Path) -> dict[str, Any] | None:
    receipt_path = output / "run_receipt.json"
    required = [output / "quant.sf", output / "cmd_info.json", output / "aux_info/meta_info.json"]
    if output.is_dir() and receipt_path.is_file() and all(path.is_file() for path in required):
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("status") == "completed":
            return receipt
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--array-index", type=int, required=True)
    parser.add_argument("--salmon", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--fastq-root", type=Path, required=True)
    parser.add_argument("--quant-root", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--library-type", default="A")
    args = parser.parse_args()
    if args.threads < 1:
        raise ValueError("threads must be positive")

    runs = load_rna_run_specs(args.manifest, study_accession=args.study)
    if args.array_index < 0 or args.array_index >= len(runs):
        raise IndexError(f"array index {args.array_index} outside 0..{len(runs) - 1}")
    run = runs[args.array_index]
    output = args.quant_root / run.study_accession / run.run_accession
    existing = _validate_existing_quant(output)
    if existing is not None:
        print(json.dumps({"action": "already_complete", "receipt": existing}, indent=2))
        return
    if output.exists():
        raise ValueError(f"Incomplete existing quantification requires inspection: {output}")
    if not args.index.is_dir():
        raise FileNotFoundError(f"Salmon index is missing: {args.index}")

    fastq_dir = args.fastq_root / run.study_accession / run.run_accession
    fastq_dir.mkdir(parents=True, exist_ok=True)
    fastq_receipts = []
    fastq_paths = []
    for spec in run.fastqs:
        target = fastq_dir / f"{run.run_accession}_{spec.mate}.fastq.gz"
        fastq_receipts.append(_download_fastq(spec, target))
        fastq_paths.append(target)

    output.parent.mkdir(parents=True, exist_ok=True)
    staging_root = Path(tempfile.mkdtemp(prefix=f".{run.run_accession}.", dir=output.parent))
    staging_output = staging_root / "quant"
    command = [
        str(args.salmon),
        "quant",
        "--index",
        str(args.index),
        "--libType",
        args.library_type,
        "--mates1",
        str(fastq_paths[0]),
        "--mates2",
        str(fastq_paths[1]),
        "--threads",
        str(args.threads),
        "--validateMappings",
        "--seqBias",
        "--gcBias",
        "--output",
        str(staging_output),
    ]
    try:
        _run(command)
        meta_path = staging_output / "aux_info/meta_info.json"
        quant_path = staging_output / "quant.sf"
        if not meta_path.is_file() or not quant_path.is_file() or quant_path.stat().st_size == 0:
            raise RuntimeError("Salmon returned without complete quantification files")
        metadata = json.loads(meta_path.read_text(encoding="utf-8"))
        receipt = {
            "status": "completed",
            "completed_at_utc": datetime.now(UTC).isoformat(),
            "repo_commit": _run(["git", "rev-parse", "HEAD"], capture=True),
            "task_index": args.array_index,
            "run": run.to_dict(),
            "fastqs": fastq_receipts,
            "salmon_version": _run([str(args.salmon), "--version"], capture=True),
            "salmon_command": command,
            "library_type_requested": args.library_type,
            "library_types_detected": metadata.get("library_types"),
            "num_processed": metadata.get("num_processed"),
            "num_mapped": metadata.get("num_mapped"),
            "percent_mapped": metadata.get("percent_mapped"),
        }
        write_json(staging_output / "run_receipt.json", receipt)
        os.replace(staging_output, output)
        staging_root.rmdir()
    except Exception as error:
        raise RuntimeError(
            f"Quantification staging directory retained: {staging_root}"
        ) from error
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
