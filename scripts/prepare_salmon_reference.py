#!/usr/bin/env python3
"""Fetch, verify, and index the frozen GENCODE v32 Salmon reference."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

from hgsoc_corneto.io import write_json
from hgsoc_corneto.rna import file_md5, resumable_curl_command


def _run(command: list[str], *, capture: bool = False) -> str:
    result = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
    )
    return result.stdout.strip() if capture else ""


def _download(url: str, target: Path, expected_md5: str) -> dict[str, Any]:
    if target.exists():
        actual_md5 = file_md5(target)
        if actual_md5 != expected_md5:
            raise ValueError(f"Existing reference has wrong MD5: {target} ({actual_md5})")
    else:
        partial = target.with_name(target.name + ".partial")
        _run(resumable_curl_command(url=url, target=partial))
        actual_md5 = file_md5(partial)
        if actual_md5 != expected_md5:
            raise ValueError(f"Downloaded reference has wrong MD5: {partial} ({actual_md5})")
        os.replace(partial, target)
    return {
        "url": url,
        "path": str(target),
        "bytes": target.stat().st_size,
        "md5": expected_md5,
    }


def _write_decoys(genome_gz: Path, target: Path) -> int:
    partial = target.with_name(target.name + ".partial")
    count = 0
    with gzip.open(genome_gz, "rt", encoding="ascii") as source, partial.open(
        "w", encoding="ascii", newline=""
    ) as output:
        for line in source:
            if line.startswith(">"):
                output.write(line[1:].split(maxsplit=1)[0] + "\n")
                count += 1
    os.replace(partial, target)
    return count


def _write_gentrome(transcripts_gz: Path, genome_gz: Path, target: Path) -> None:
    partial = target.with_name(target.name + ".partial")
    with partial.open("wb") as output:
        for source_path in (transcripts_gz, genome_gz):
            with gzip.open(source_path, "rb") as source:
                shutil.copyfileobj(source, output, length=16 * 1024 * 1024)
    os.replace(partial, target)


def _index_is_complete(index_dir: Path) -> bool:
    required = ("versionInfo.json", "info.json", "seq.bin")
    return index_dir.is_dir() and all((index_dir / name).is_file() for name in required)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--salmon", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.threads < 1:
        raise ValueError("threads must be positive")

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    reference_config = config["reference"]
    reference_root = args.reference_root or Path(config["roihu_paths"]["reference"])
    downloads = reference_root / "downloads"
    downloads.mkdir(parents=True, exist_ok=True)
    base_url = reference_config["base_url"].rstrip("/")

    fetched = {}
    paths = {}
    for key, source in reference_config["files"].items():
        target = downloads / source["name"]
        fetched[key] = _download(f"{base_url}/{source['name']}", target, source["md5"])
        paths[key] = target

    decoys = reference_root / "decoys.txt"
    gentrome = reference_root / "gentrome.fa"
    decoy_count = _write_decoys(paths["genome"], decoys) if not decoys.exists() else sum(
        1 for _ in decoys.open(encoding="ascii")
    )
    if not gentrome.exists():
        _write_gentrome(paths["transcripts"], paths["genome"], gentrome)

    index_dir = reference_root / "salmon_index"
    command = [
        str(args.salmon),
        "index",
        "--transcripts",
        str(gentrome),
        "--decoys",
        str(decoys),
        "--index",
        str(index_dir),
        "--threads",
        str(args.threads),
        "--kmerLen",
        str(reference_config["index"]["kmer_length"]),
        "--gencode",
    ]
    if not _index_is_complete(index_dir):
        if index_dir.exists():
            raise ValueError(f"Incomplete existing index requires inspection: {index_dir}")
        staging_root = Path(tempfile.mkdtemp(prefix=".salmon-index.", dir=reference_root))
        staging_index = staging_root / "index"
        staging_command = list(command)
        staging_command[staging_command.index(str(index_dir))] = str(staging_index)
        try:
            _run(staging_command)
            if not _index_is_complete(staging_index):
                raise RuntimeError("Salmon returned without a complete index")
            os.replace(staging_index, index_dir)
            staging_root.rmdir()
        except Exception as error:
            raise RuntimeError(
                f"Index staging directory retained for audit: {staging_root}"
            ) from error

    version_info = json.loads((index_dir / "versionInfo.json").read_text(encoding="utf-8"))
    receipt = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "repo_commit": _run(["git", "rev-parse", "HEAD"], capture=True),
        "salmon_version": _run([str(args.salmon), "--version"], capture=True),
        "reference": fetched,
        "gentrome": {"path": str(gentrome), "bytes": gentrome.stat().st_size},
        "decoys": {"path": str(decoys), "count": decoy_count},
        "index": {"path": str(index_dir), "version_info": version_info},
        "index_command": command,
    }
    write_json(args.receipt, receipt)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
