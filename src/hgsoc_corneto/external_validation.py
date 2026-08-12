"""Fail-closed contracts for independent organoid and DepMap validation.

The helpers in this module deliberately separate *model selection* from matrix
extraction.  In particular, an ovarian lineage label is never treated as proof
that a DepMap model faithfully represents HGSOC.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import os
import re
import shutil
import statistics
import tempfile
import urllib.request
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any, TextIO

_GENE_LABEL = re.compile(r"^(?P<symbol>.+?)\s+\(\d+\)$")


def _open_text(path: str | Path) -> TextIO:
    target = Path(path)
    if target.suffix == ".gz":
        return gzip.open(target, "rt", encoding="utf-8-sig", newline="")
    return target.open(encoding="utf-8-sig", newline="")


def file_sha256(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def download_verified(
    url: str,
    target: str | Path,
    *,
    expected_sha256: str,
    timeout_seconds: int = 180,
) -> dict[str, Any]:
    """Atomically download one immutable public file and verify it before publish."""

    destination = Path(target)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and file_sha256(destination) == expected_sha256:
        return {
            "path": str(destination),
            "sha256": expected_sha256,
            "bytes": destination.stat().st_size,
            "downloaded": False,
        }
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
        request = urllib.request.Request(url, headers={"User-Agent": "hgsoc-corneto/0.1"})
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                shutil.copyfileobj(response, handle)
            # NamedTemporaryFile remains open until this context exits.  Hashing
            # its path before explicitly flushing can observe an empty or
            # truncated file on some parallel filesystems (as seen on Roihu).
            handle.flush()
            os.fsync(handle.fileno())
            observed = file_sha256(temporary)
            if observed != expected_sha256:
                raise ValueError(
                    f"download checksum mismatch: expected {expected_sha256}, observed {observed}"
                )
            temporary.replace(destination)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    return {
        "path": str(destination),
        "sha256": expected_sha256,
        "bytes": destination.stat().st_size,
        "downloaded": True,
    }


def depmap_download_preflight(
    *,
    release: str | None,
    model_path: str | Path | None,
    gene_effect_path: str | Path | None,
    release_readme_path: str | Path | None,
    landing_url: str,
) -> dict[str, Any]:
    """Return an auditable ready/blocked state without inventing download URLs."""

    reasons: list[str] = []
    release_pattern = re.compile(r"^\d{2}Q[1-4]$")
    if release is None or not release_pattern.fullmatch(release):
        reasons.append("an explicit quarterly DepMap release such as 26Q1 is required")

    declared = {
        "Model.csv": Path(model_path) if model_path is not None else None,
        "CRISPRGeneEffect.csv": (
            Path(gene_effect_path) if gene_effect_path is not None else None
        ),
        "release_README": (
            Path(release_readme_path) if release_readme_path is not None else None
        ),
    }
    files: dict[str, dict[str, Any]] = {}
    for label, path in declared.items():
        if path is None:
            reasons.append(f"{label} was not supplied")
        elif not path.is_file():
            reasons.append(f"{label} does not exist at {path}")
        else:
            files[label] = {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }

    blocked = bool(reasons)
    return {
        "scientific_success": False,
        "status": "blocked" if blocked else "ready_for_schema_validation",
        "release": release,
        "official_download_landing_page": landing_url,
        "files": files,
        "blocking_reasons": reasons,
        "claim_limit": (
            "This preflight is not a dependency result. No direct data URL is guessed, and files "
            "from unspecified or mixed DepMap releases are not accepted."
        ),
    }


def audit_gse_count_matrix(
    path: str | Path,
    *,
    expected_sha256: str,
    expected_samples: Iterable[str],
    expected_gene_rows: int,
) -> dict[str, Any]:
    """Validate the immutable GEO count-matrix payload and its tabular schema."""

    observed_sha256 = file_sha256(path)
    if observed_sha256 != expected_sha256:
        raise ValueError(
            f"GEO checksum mismatch: expected {expected_sha256}, observed {observed_sha256}"
        )

    expected = list(expected_samples)
    genes: set[str] = set()
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("GEO count matrix is empty") from exc
        if not header or header[0] != "":
            raise ValueError("GEO count matrix must have an empty gene-ID header cell")
        samples = header[1:]
        if set(samples) != set(expected) or len(samples) != len(expected):
            raise ValueError(
                f"GEO sample set mismatch: expected={sorted(expected)}; observed={sorted(samples)}"
            )
        if len(samples) != len(set(samples)):
            raise ValueError("GEO count matrix contains duplicate sample columns")

        row_count = 0
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"GEO row {line_number} has {len(row)} fields; expected {len(header)}"
                )
            gene_id = row[0]
            if not gene_id or gene_id in genes:
                raise ValueError(f"GEO row {line_number} has a missing or duplicate gene ID")
            genes.add(gene_id)
            for value in row[1:]:
                try:
                    count = int(value)
                except ValueError as exc:
                    raise ValueError(f"GEO row {line_number} contains a non-integer count") from exc
                if count < 0:
                    raise ValueError(f"GEO row {line_number} contains a negative count")
            row_count += 1

    if row_count != expected_gene_rows:
        raise ValueError(
            f"GEO gene-row mismatch: expected {expected_gene_rows}, observed {row_count}"
        )
    return {
        "scientific_success": True,
        "sha256": observed_sha256,
        "sample_count": len(samples),
        "gene_row_count": row_count,
        "samples": samples,
    }


def extract_gse_candidate_log_cpm(
    matrix_path: str | Path,
    *,
    sample_groups: dict[str, str],
    gene_id_to_symbol: dict[str, str],
    candidate_symbols: Iterable[str],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Extract candidate expression without inventing an Ensembl-symbol mapping.

    Counts for multiple Ensembl IDs explicitly mapped to the same symbol are
    summed before library-size normalization.  The returned values are
    ``log2(CPM + 1)`` and are intended for descriptive external validation.
    """

    candidates = {symbol.strip().upper() for symbol in candidate_symbols if symbol.strip()}
    if not candidates:
        raise ValueError("at least one candidate gene symbol is required")
    normalized_map: dict[str, str] = {}
    for gene_id, symbol in gene_id_to_symbol.items():
        normalized_id = gene_id.strip().split(".", maxsplit=1)[0]
        normalized_symbol = symbol.strip().upper()
        if not normalized_id or not normalized_symbol:
            raise ValueError("gene mapping contains an empty gene ID or symbol")
        previous = normalized_map.get(normalized_id)
        if previous is not None and previous != normalized_symbol:
            raise ValueError(f"conflicting symbols for Ensembl gene ID {normalized_id}")
        normalized_map[normalized_id] = normalized_symbol

    with _open_text(matrix_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        if set(samples) != set(sample_groups):
            raise ValueError("GSE sample groups do not match matrix columns")
        library_sizes = [0] * len(samples)
        candidate_counts = {symbol: [0] * len(samples) for symbol in candidates}
        observed_symbols: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"GSE row {line_number} has the wrong number of fields")
            counts = [int(value) for value in row[1:]]
            for index, count in enumerate(counts):
                if count < 0:
                    raise ValueError(f"GSE row {line_number} contains a negative count")
                library_sizes[index] += count
            symbol = normalized_map.get(row[0].split(".", maxsplit=1)[0])
            if symbol in candidates:
                observed_symbols.add(symbol)
                for index, count in enumerate(counts):
                    candidate_counts[symbol][index] += count

    if any(size <= 0 for size in library_sizes):
        raise ValueError("GSE matrix contains an empty sample library")
    output = []
    for symbol in sorted(observed_symbols):
        for index, sample in enumerate(samples):
            count = candidate_counts[symbol][index]
            log_cpm = math.log2((count * 1_000_000 / library_sizes[index]) + 1.0)
            output.append(
                {
                    "gene_symbol": symbol,
                    "sample_id": sample,
                    "validation_group": sample_groups[sample],
                    "raw_count": count,
                    "library_size": library_sizes[index],
                    "log2_cpm_plus_1": log_cpm,
                }
            )
    return output, sorted(candidates - observed_symbols)


def summarize_expression_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Summarize model-level log2(CPM+1) values without pseudo-replication tests."""

    buckets: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        key = (str(row["gene_symbol"]), str(row["validation_group"]))
        buckets.setdefault(key, []).append(float(row["log2_cpm_plus_1"]))
    return [
        {
            "gene_symbol": gene_symbol,
            "validation_group": group,
            "n_models": len(values),
            "mean_log2_cpm_plus_1": statistics.fmean(values),
            "median_log2_cpm_plus_1": statistics.median(values),
            "min_log2_cpm_plus_1": min(values),
            "max_log2_cpm_plus_1": max(values),
        }
        for (gene_symbol, group), values in sorted(buckets.items())
    ]


def normalize_cell_line_name(value: str) -> str:
    """Normalize punctuation and optional CCLE lineage suffixes for exact matching."""

    base = value.strip().upper()
    for suffix in ("_OVARY", "_FALLOPIAN_TUBE"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return re.sub(r"[^A-Z0-9]", "", base)


def parse_depmap_gene_label(value: str) -> str:
    """Return an HGNC-like symbol from ``SYMBOL (EntrezID)`` or ``SYMBOL``."""

    match = _GENE_LABEL.fullmatch(value.strip())
    return (match.group("symbol") if match else value).strip().upper()


def read_tsv_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv_rows(path: str | Path) -> list[dict[str, str]]:
    with _open_text(path) as handle:
        return list(csv.DictReader(handle))


def resolve_depmap_hgsoc_models(
    model_rows: list[dict[str, str]], curated_rows: list[dict[str, str]]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Resolve a literature-curated positive set and an ovarian comparator set.

    Model IDs are release-specific assertions.  A mismatch between the expected
    ID and model name fails closed instead of silently selecting a different
    line.  Other ovarian models are returned only as comparators and are never
    relabelled as HGSOC.
    """

    required_model_columns = {"ModelID", "StrippedCellLineName", "OncotreeLineage"}
    if not model_rows:
        raise ValueError("DepMap Model.csv contains no rows")
    missing_columns = required_model_columns - set(model_rows[0])
    if missing_columns:
        raise ValueError(f"DepMap Model.csv missing columns: {sorted(missing_columns)}")
    required_curated = {"model_id", "cell_line", "evidence_tier", "evidence_url"}
    if not curated_rows or required_curated - set(curated_rows[0]):
        raise ValueError("curated model table is empty or missing required columns")

    by_id = {row["ModelID"]: row for row in model_rows}
    selected: list[dict[str, str]] = []
    selected_ids: set[str] = set()
    for assertion in curated_rows:
        model_id = assertion["model_id"]
        if model_id in selected_ids:
            raise ValueError(f"duplicate curated DepMap model ID: {model_id}")
        model = by_id.get(model_id)
        if model is None:
            raise ValueError(f"curated DepMap model absent from this release: {model_id}")
        observed = normalize_cell_line_name(model["StrippedCellLineName"])
        expected = normalize_cell_line_name(assertion["cell_line"])
        if observed != expected:
            raise ValueError(
                f"DepMap identity mismatch for {model_id}: expected {expected}, observed {observed}"
            )
        if assertion["evidence_tier"] != "high_confidence_hgsoc_like":
            raise ValueError(f"unsupported positive-set evidence tier for {model_id}")
        selected_ids.add(model_id)
        selected.append({**model, **assertion, "validation_group": "hgsoc_like_positive"})

    ovarian_comparators = [
        {**row, "validation_group": "other_ovarian_not_hgsoc_positive"}
        for row in model_rows
        if row["ModelID"] not in selected_ids
        and normalize_cell_line_name(row.get("OncotreeLineage", ""))
        in {"OVARYFALLOPIANTUBE", "OVARY"}
    ]
    return selected, ovarian_comparators


def extract_depmap_gene_effects(
    matrix_path: str | Path,
    *,
    model_groups: dict[str, str],
    candidate_symbols: Iterable[str],
) -> tuple[list[dict[str, Any]], list[str], list[str]]:
    """Stream a DepMap gene-effect matrix and return a tidy selected subset."""

    candidates = {symbol.strip().upper() for symbol in candidate_symbols if symbol.strip()}
    if not candidates:
        raise ValueError("at least one candidate gene symbol is required")

    observed_models: set[str] = set()
    with _open_text(matrix_path) as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("DepMap gene-effect matrix is empty") from exc
        if len(header) < 2:
            raise ValueError("DepMap gene-effect matrix has no gene columns")
        symbols = [parse_depmap_gene_label(label) for label in header[1:]]
        if len(symbols) != len(set(symbols)):
            raise ValueError("DepMap gene-effect matrix has duplicate normalized gene symbols")
        selected_indexes = [index for index, symbol in enumerate(symbols) if symbol in candidates]
        observed_candidates = {symbols[index] for index in selected_indexes}

        output: list[dict[str, Any]] = []
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"DepMap gene-effect row {line_number} has {len(row)} fields; "
                    f"expected {len(header)}"
                )
            model_id = row[0]
            group = model_groups.get(model_id)
            if group is None:
                continue
            observed_models.add(model_id)
            for index in selected_indexes:
                raw_value = row[index + 1]
                if raw_value == "":
                    continue
                value = float(raw_value)
                if not math.isfinite(value):
                    raise ValueError(
                        f"non-finite DepMap value for {model_id}, {symbols[index]}"
                    )
                output.append(
                    {
                        "model_id": model_id,
                        "validation_group": group,
                        "gene_symbol": symbols[index],
                        "gene_effect": value,
                    }
                )

    missing_models = sorted(set(model_groups) - observed_models)
    missing_genes = sorted(candidates - observed_candidates)
    return output, missing_models, missing_genes


def summarize_dependency_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Produce descriptive, non-causal summaries by gene and declared group."""

    buckets: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        key = (str(row["gene_symbol"]), str(row["validation_group"]))
        buckets.setdefault(key, []).append(float(row["gene_effect"]))
    summaries = []
    for (gene_symbol, group), values in sorted(buckets.items()):
        summaries.append(
            {
                "gene_symbol": gene_symbol,
                "validation_group": group,
                "n_models": len(values),
                "mean_gene_effect": statistics.fmean(values),
                "median_gene_effect": statistics.median(values),
                "min_gene_effect": min(values),
                "max_gene_effect": max(values),
            }
        )
    return summaries


def iter_candidate_symbols(path: str | Path) -> Iterator[str]:
    rows = read_tsv_rows(path)
    if not rows or "gene_symbol" not in rows[0]:
        raise ValueError("candidate TSV must contain a gene_symbol column")
    for row in rows:
        symbol = row["gene_symbol"].strip()
        if symbol:
            yield symbol


def read_gene_map(path: str | Path) -> dict[str, str]:
    rows = read_tsv_rows(path)
    required = {"gene_id", "gene_symbol"}
    if not rows or required - set(rows[0]):
        raise ValueError("gene-map TSV must contain gene_id and gene_symbol columns")
    mapping: dict[str, str] = {}
    for row in rows:
        gene_id = row["gene_id"].strip()
        symbol = row["gene_symbol"].strip()
        if not gene_id or not symbol:
            raise ValueError("gene-map TSV contains an empty gene_id or gene_symbol")
        previous = mapping.get(gene_id)
        if previous is not None and previous != symbol:
            raise ValueError(f"conflicting gene-map rows for {gene_id}")
        mapping[gene_id] = symbol
    return mapping
