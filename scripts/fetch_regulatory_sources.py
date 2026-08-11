#!/usr/bin/env python3
"""Fetch small, public regulatory priors and write provenance receipts.

The downloaded tables are intentionally kept outside Git.  This command
records the exact URL, retrieval time, byte count, SHA-256, declared licence
note, and tabular column statistics.  It does not query a solver or read any
phenotype/licence file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError("PyYAML is required to read the regulatory source config") from error
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a mapping in {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _column_stats(path: Path, required_any: list[list[str]]) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"Empty regulatory table: {path}") from error
        if not header or len(header) != len(set(header)):
            raise ValueError(f"Missing or duplicate header columns: {path}")
        rows = sum(1 for _ in reader)
    missing_groups = [group for group in required_any if not set(group).intersection(header)]
    return {
        "rows": rows,
        "columns": len(header),
        "column_names": header,
        "missing_required_any_groups": missing_groups,
    }


def _download(
    *,
    url: str,
    destination: Path,
    user_agent: str,
    timeout: int,
    max_bytes: int,
    overwrite: bool,
) -> dict[str, Any]:
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing source: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".partial")
    if partial.exists():
        partial.unlink()
    started = datetime.now(timezone.utc)
    request = urllib.request.Request(url, headers={"User-Agent": user_agent})
    byte_count = 0
    with urllib.request.urlopen(request, timeout=timeout) as response, partial.open("wb") as handle:
        content_length = response.headers.get("Content-Length")
        declared_bytes = int(content_length) if content_length and content_length.isdigit() else None
        if declared_bytes is not None and declared_bytes > max_bytes:
            raise ValueError(f"Refusing {url}: Content-Length {declared_bytes} exceeds {max_bytes}")
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            byte_count += len(chunk)
            if byte_count > max_bytes:
                raise ValueError(f"Refusing {url}: download exceeds {max_bytes} bytes")
            handle.write(chunk)
    os.replace(partial, destination)
    finished = datetime.now(timezone.utc)
    return {
        "url": url,
        "retrieved_at_utc": finished.isoformat(),
        "started_at_utc": started.isoformat(),
        "bytes": byte_count,
        "sha256": _sha256(destination),
        "content_type": response.headers.get("Content-Type"),
        "declared_content_length": declared_bytes,
    }


def fetch(
    *, config_path: Path, output_dir: Path, overwrite: bool = False
) -> dict[str, Any]:
    config = _load_yaml(config_path)
    retrieval = config.get("retrieval", {})
    sources = config.get("sources")
    if not isinstance(sources, dict) or not sources:
        raise ValueError("Config must contain a non-empty sources mapping")
    output_dir.mkdir(parents=True, exist_ok=True)
    receipt: dict[str, Any] = {
        "schema_version": config.get("schema_version"),
        "status": "completed",
        "config_path": str(config_path),
        "sources": {},
        "response_blind": True,
    }
    failures: list[str] = []
    for name, raw in sources.items():
        if not isinstance(raw, dict):
            failures.append(f"{name}: source config is not a mapping")
            continue
        destination = output_dir / str(raw.get("output", f"{name}.tsv"))
        record: dict[str, Any] = {
            "role": raw.get("role"),
            "license": raw.get("license"),
            "path": str(destination),
            "url": raw.get("url"),
        }
        try:
            url = str(raw["url"])
            transfer = _download(
                url=url,
                destination=destination,
                user_agent=str(retrieval.get("user_agent", "hgsoc-corneto-regulatory")),
                timeout=int(retrieval.get("timeout_seconds", 120)),
                max_bytes=int(retrieval.get("max_bytes", 50_000_000)),
                overwrite=overwrite,
            )
            required_any = [list(group) for group in raw.get("required_any_columns", [])]
            columns = _column_stats(destination, required_any)
            if columns["missing_required_any_groups"]:
                raise ValueError(
                    f"{name}: required column groups missing: "
                    f"{columns['missing_required_any_groups']}"
                )
            record.update(transfer)
            record["column_stats"] = columns
            record["status"] = "downloaded"
        except Exception as error:  # preserve a machine-readable blocked receipt
            record.update({"status": "blocked", "error": f"{type(error).__name__}: {error}"})
            failures.append(f"{name}: {error}")
        receipt["sources"][name] = record

    if failures:
        receipt["status"] = "blocked_download"
        receipt["errors"] = failures
    receipt_path = output_dir / "regulatory_sources_receipt.json"
    if receipt_path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite receipt: {receipt_path}")
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=output_dir, delete=False
    ) as handle:
        temporary = Path(handle.name)
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, receipt_path)
    if failures:
        raise RuntimeError("; ".join(failures))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/regulatory_sources.yaml"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    receipt = fetch(config_path=args.config, output_dir=args.output_dir, overwrite=args.overwrite)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
