#!/usr/bin/env python3
"""Fetch and verify the pinned public inputs for the Meeson benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path

MEESON_COMMIT = "1b0e9050720be8fb6769a9a7ab14a4b19bbce8c3"
MEESON_URL = f"https://codeload.github.com/katemeeson/PhD_2024/tar.gz/{MEESON_COMMIT}"
MEESON_SHA256 = "dbd5a69aac4a5d1a0f35298cd2f87494edad54f561b4a238c440a418aecac5a0"

HUMAN_GEM_RELEASE = "v1.4.1"
HUMAN_GEM_URL = (
    "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/"
    "v1.4.1/model/Human-GEM.xml"
)
HUMAN_GEM_SHA256 = "57d1b137f0c90d83a3e4f9a8225d74d37523594e6ee99f622b160a014d9f7050"


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verified_download(url: str, destination: Path, expected_sha256: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and sha256(destination) == expected_sha256:
        return
    partial = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "hgsoc-corneto/0.1"})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    observed = sha256(partial)
    if observed != expected_sha256:
        partial.unlink(missing_ok=True)
        raise ValueError(
            f"Checksum mismatch for {url}: expected {expected_sha256}, observed {observed}"
        )
    partial.replace(destination)


def safe_extract_tar(archive: Path, destination: Path) -> Path:
    expected_root = f"PhD_2024-{MEESON_COMMIT}"
    extracted_root = destination / expected_root
    marker = extracted_root / "src" / "integrate_omics.py"
    if marker.is_file():
        return extracted_root

    destination.mkdir(parents=True, exist_ok=True)
    resolved_destination = destination.resolve()
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            target = (destination / member.name).resolve()
            if not target.is_relative_to(resolved_destination):
                raise ValueError(f"Unsafe archive member: {member.name}")
        bundle.extractall(destination, filter="data")
    if not marker.is_file():
        raise FileNotFoundError(f"Expected file absent after extraction: {marker}")
    return extracted_root


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()

    destination = args.destination.resolve()
    downloads = destination / "downloads"
    human_gem = downloads / "Human-GEM-v1.4.1.xml"
    meeson_archive = downloads / f"PhD_2024-{MEESON_COMMIT}.tar.gz"
    verified_download(HUMAN_GEM_URL, human_gem, HUMAN_GEM_SHA256)
    verified_download(MEESON_URL, meeson_archive, MEESON_SHA256)
    checkout = safe_extract_tar(meeson_archive, destination / "sources")

    receipt = {
        "human_gem": {
            "release": HUMAN_GEM_RELEASE,
            "url": HUMAN_GEM_URL,
            "path": str(human_gem),
            "bytes": human_gem.stat().st_size,
            "sha256": sha256(human_gem),
        },
        "meeson_phd_2024": {
            "commit": MEESON_COMMIT,
            "url": MEESON_URL,
            "archive_path": str(meeson_archive),
            "archive_bytes": meeson_archive.stat().st_size,
            "archive_sha256": sha256(meeson_archive),
            "checkout": str(checkout),
        },
    }
    (destination / "source_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
