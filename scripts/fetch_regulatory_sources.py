#!/usr/bin/env python3
"""Fetch and normalize public CollecTRI and OmniPath interaction tables.

The raw and normalized tables are intentionally external to Git.  This script
creates a provenance receipt containing URLs, retrieval time, byte counts,
SHA256 digests, source columns and filtering counts.  It does not read any
phenotype or license-file content.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_COLLECTRI = (
    "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&format=tsv"
)
DEFAULT_OMNIPATH = (
    "https://omnipathdb.org/interactions?datasets=omnipath&genesymbols=yes&format=tsv"
)
REQUIRED = (
    "source_genesymbol",
    "target_genesymbol",
    "is_directed",
    "consensus_stimulation",
    "consensus_inhibition",
)


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _sha256(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _download(url: str, destination: Path) -> dict[str, object]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "hgsoc-corneto-regulatory/1.0", "Accept-Encoding": "identity"},
    )
    part = destination.with_suffix(destination.suffix + ".part")
    digest = hashlib.sha256()
    size = 0
    with urllib.request.urlopen(request, timeout=120) as response, part.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            digest.update(chunk)
            size += len(chunk)
    os.replace(part, destination)
    return {
        "url": url,
        "http_status": getattr(response, "status", None),
        "content_type": response.headers.get("content-type"),
        "etag": response.headers.get("etag"),
        "last_modified": response.headers.get("last-modified"),
        "bytes": size,
        "sha256": digest.hexdigest(),
    }


def _normalize(raw: Path, output: Path, source_name: str) -> dict[str, object]:
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = {"input_rows": 0, "kept_rows": 0, "skipped_missing": 0,
              "skipped_undirected": 0, "skipped_ambiguous_sign": 0,
              "skipped_self_loop": 0}
    unique: set[tuple[str, str, int]] = set()
    columns: list[str] = []
    with raw.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = list(reader.fieldnames or [])
        missing = [column for column in REQUIRED if column not in columns]
        if missing:
            raise ValueError(f"{source_name}: missing required columns: {missing}")
        for row in reader:
            counts["input_rows"] += 1
            source = (row.get("source_genesymbol") or "").strip()
            target = (row.get("target_genesymbol") or "").strip()
            if not source or not target:
                counts["skipped_missing"] += 1
                continue
            if not _truthy(row.get("is_directed")):
                counts["skipped_undirected"] += 1
                continue
            stimulation = _truthy(row.get("consensus_stimulation"))
            inhibition = _truthy(row.get("consensus_inhibition"))
            if stimulation == inhibition:
                counts["skipped_ambiguous_sign"] += 1
                continue
            if source == target:
                counts["skipped_self_loop"] += 1
                continue
            unique.add((source, target, 1 if stimulation else -1))
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source", "target", "sign"])
        for source, target, sign in sorted(unique):
            writer.writerow([source, target, sign])
    counts["kept_rows"] = len(unique)
    digest, size = _sha256(output)
    return {"source": source_name, "columns": columns, **counts,
            "normalized_path": str(output), "normalized_bytes": size,
            "normalized_sha256": digest}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--collectri-url", default=DEFAULT_COLLECTRI)
    parser.add_argument("--omnipath-url", default=DEFAULT_OMNIPATH)
    parser.add_argument("--keep-raw", action="store_true")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    retrieved = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    receipt: dict[str, object] = {
        "schema_version": "regulatory_sources.v1",
        "status": "running",
        "organism": "human",
        "retrieved_at_utc": retrieved,
        "normalization": {
            "require_directed": True,
            "require_unambiguous_consensus_sign": True,
            "remove_self_loops": True,
        },
        "sources": {},
    }
    raw_paths: list[Path] = []
    try:
        for name, url, normalized_name in (
            ("collectri", args.collectri_url, "collectri_tf_target.tsv.gz"),
            ("omnipath", args.omnipath_url, "omnipath_signed_pkn.tsv.gz"),
        ):
            raw = args.out_dir / f"{name}.raw.tsv"
            raw_paths.append(raw)
            download = _download(url, raw)
            normalized = _normalize(raw, args.out_dir / normalized_name, name)
            receipt["sources"][name] = {**download, **normalized,
                "license_note": (
                    "CollecTRI GPL-3 and original resource terms apply."
                    if name == "collectri" else
                    "OmniPath resource-specific licensing and attribution terms apply."
                )}
        receipt["status"] = "completed"
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error_class"] = type(exc).__name__
        receipt["error_message"] = str(exc)[:500]
        raise
    finally:
        receipt_path = args.out_dir / "regulatory_source_receipt.json"
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if not args.keep_raw:
            for raw in raw_paths:
                raw.unlink(missing_ok=True)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
