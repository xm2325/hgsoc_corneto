"""Deterministic tabular and provenance I/O."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any, TextIO


def _serialize(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def write_tsv(
    path: str | Path, rows: list[dict[str, Any]], fields: list[str] | None = None
) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not rows and fields is None:
        raise ValueError("fields are required for an empty table")
    columns = fields or list(rows[0])
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _serialize(row.get(column)) for column in columns})


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def checksum_rows(paths: Iterable[Path], root: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(paths):
        rows.append(
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return rows


def write_json(path: str | Path, value: Any) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


@contextmanager
def deterministic_gzip_text_writer(path: str | Path) -> Iterator[TextIO]:
    """Write gzip text without timestamps or filenames in the gzip header."""

    target = Path(path)
    with target.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_handle:
                yield text_handle
