"""Provenance-checked preparation of the paired-site GSE277107 RNA dataset.

The GEO series describes matched ovary/adnexal and omentum HGSC tissue from
11 people.  GEO does not expose a dedicated ``patient_id`` characteristic.
Consequently the code below derives a *pair key* only from the shared numeric
prefix in the submitted sample descriptions, and records that derivation in
every output row.  It must not be treated as a clinical identifier.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import shutil
import tempfile
import urllib.request
from collections import defaultdict
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TextIO

from hgsoc_corneto.io import deterministic_gzip_text_writer, sha256, write_json, write_tsv


class GSE277107Error(ValueError):
    """Raised when public metadata or expression inputs fail a frozen gate."""


PAIR_KEY_DERIVATION = (
    "shared WM+numeric prefix parsed from GEO Sample_description; pairing is supported by "
    "GSE277107 Series_overall_design, not by a standalone patient_id field"
)

_DESCRIPTION = re.compile(r"^WM(?P<pair>[0-9]+)(?P<aliquot>[A-Z])_(?P<site>OV|OM)$")
_TITLE = re.compile(
    r"^(?P<pair>[0-9]+)(?P<aliquot>[A-Z])_(?P<site>OV|OM)_"
    r"(?P<histology>HGSC|HGSCS)_(?P<tissue>ovary|omentum)$"
)


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def load_source_manifest(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != "hgsoc_external_sources.v1":
        raise GSE277107Error("unsupported source-manifest schema")
    if data.get("study_accession") != "GSE277107":
        raise GSE277107Error("source manifest is not for GSE277107")
    files = data.get("files")
    if (
        not isinstance(files, list)
        or len(files) != 2
        or {item.get("role") for item in files}
        != {"geo_family_soft", "processed_tpm"}
    ):
        raise GSE277107Error("source manifest must define exactly SOFT and TPM roles")
    for item in files:
        checksum = str(item.get("sha256", ""))
        if not re.fullmatch(r"[0-9a-f]{64}", checksum):
            raise GSE277107Error(f"invalid SHA256 for {item.get('role')}")
        if not str(item.get("url", "")).startswith("https://ftp.ncbi.nlm.nih.gov/geo/"):
            raise GSE277107Error(f"non-NCBI source URL for {item.get('role')}")
        if Path(str(item.get("filename", ""))).name != item.get("filename"):
            raise GSE277107Error(f"unsafe filename for {item.get('role')}")
    return data


def verify_source_files(
    source_manifest: dict[str, Any], source_dir: Path
) -> dict[str, Path]:
    resolved: dict[str, Path] = {}
    for item in source_manifest["files"]:
        path = source_dir / item["filename"]
        if not path.is_file():
            raise GSE277107Error(f"missing {item['role']} source file: {path}")
        observed_bytes = path.stat().st_size
        observed_sha = sha256(path)
        if observed_bytes != int(item["bytes"]):
            raise GSE277107Error(
                f"byte-size mismatch for {item['role']}: {observed_bytes} != {item['bytes']}"
            )
        if observed_sha != item["sha256"]:
            raise GSE277107Error(f"SHA256 mismatch for {item['role']}: {observed_sha}")
        resolved[item["role"]] = path
    return resolved


def fetch_sources(source_manifest_path: Path, output_dir: Path) -> dict[str, Path]:
    """Download the two frozen GEO inputs and verify bytes plus SHA256."""

    manifest = load_source_manifest(source_manifest_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    for item in manifest["files"]:
        target = output_dir / item["filename"]
        if target.is_file() and target.stat().st_size == int(item["bytes"]):
            if sha256(target) == item["sha256"]:
                continue
        with tempfile.NamedTemporaryFile(dir=output_dir, delete=False) as handle:
            temporary = Path(handle.name)
            request = urllib.request.Request(
                item["url"], headers={"User-Agent": "hgsoc-corneto/0.1"}
            )
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    shutil.copyfileobj(response, handle)
            except Exception:
                temporary.unlink(missing_ok=True)
                raise
        observed_sha = sha256(temporary)
        if temporary.stat().st_size != int(item["bytes"]) or observed_sha != item["sha256"]:
            temporary.unlink(missing_ok=True)
            raise GSE277107Error(
                f"download verification failed for {item['role']}: SHA256={observed_sha}"
            )
        temporary.replace(target)
    return verify_source_files(manifest, output_dir)


def _parse_soft(path: Path) -> tuple[dict[str, list[str]], list[dict[str, list[str]]]]:
    series: dict[str, list[str]] = defaultdict(list)
    samples: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("^SERIES = "):
                series["^SERIES"].append(line.split(" = ", 1)[1])
                current = None
            elif line.startswith("^SAMPLE = "):
                current = defaultdict(list)
                current["^SAMPLE"].append(line.split(" = ", 1)[1])
                samples.append(current)
            elif line.startswith("!") and " = " in line:
                key, value = line.split(" = ", 1)
                (current if current is not None else series)[key].append(value)
    return dict(series), [dict(item) for item in samples]


def _one(record: dict[str, list[str]], key: str, context: str) -> str:
    values = record.get(key, [])
    if len(values) != 1 or not values[0].strip():
        raise GSE277107Error(f"{context}: expected exactly one non-empty {key}")
    return values[0].strip()


def build_paired_metadata(
    soft_path: Path, *, expected_samples: int = 22, expected_pairs: int = 11
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Parse GEO metadata and require complete, internally consistent site pairs."""

    series, sample_records = _parse_soft(soft_path)
    accession = _one(series, "!Series_geo_accession", "series")
    if accession != "GSE277107":
        raise GSE277107Error(f"unexpected GEO accession: {accession}")
    design = _one(series, "!Series_overall_design", "series")
    design_lower = design.casefold()
    if "matched primary site" not in design_lower or "omentum" not in design_lower:
        raise GSE277107Error("series design does not substantiate matched ovary/omentum samples")
    if len(sample_records) != expected_samples:
        raise GSE277107Error(
            f"expected {expected_samples} GEO samples, observed {len(sample_records)}"
        )

    rows: list[dict[str, str]] = []
    anomalies: list[dict[str, str]] = []
    for record in sample_records:
        gsm = _one(record, "^SAMPLE", "sample")
        title = _one(record, "!Sample_title", gsm)
        description = _one(record, "!Sample_description", gsm)
        source_name = _one(record, "!Sample_source_name_ch1", gsm)
        characteristic = _one(record, "!Sample_characteristics_ch1", gsm)
        title_match = _TITLE.fullmatch(title)
        description_match = _DESCRIPTION.fullmatch(description)
        if title_match is None or description_match is None:
            raise GSE277107Error(f"{gsm}: unparseable submitted title/description")
        for field in ("pair", "aliquot", "site"):
            if title_match[field] != description_match[field]:
                raise GSE277107Error(f"{gsm}: title/description disagree on {field}")
        site_code = description_match["site"]
        reported_site = "ovary" if site_code == "OV" else "omentum"
        expected_source = f"High Grade Serous Ovarian Cancer_{reported_site}"
        if source_name != expected_source or characteristic != f"tissue: {expected_source}":
            raise GSE277107Error(f"{gsm}: submitted tissue fields disagree with {site_code}")
        if title_match["tissue"] != reported_site:
            raise GSE277107Error(f"{gsm}: title tissue disagrees with description")
        relations = record.get("!Sample_relation", [])
        sra_matches = [
            re.search(r"term=(SRX[0-9]+)$", value).group(1)
            for value in relations
            if re.search(r"term=(SRX[0-9]+)$", value)
        ]
        if len(sra_matches) != 1:
            raise GSE277107Error(f"{gsm}: expected exactly one SRA experiment relation")
        anomaly = ""
        if title_match["histology"] == "HGSCS":
            anomaly = "submitted_title_histology_token_is_HGSCS_not_HGSC"
            anomalies.append(
                {"geo_sample_accession": gsm, "field": "Sample_title", "value": title}
            )
        pair_id = f"WM{description_match['pair']}"
        rows.append(
            {
                "study_accession": accession,
                "geo_sample_accession": gsm,
                "sample_id": description,
                "pair_id": pair_id,
                "aliquot_key": description_match["aliquot"],
                "reported_site_code": site_code,
                "reported_site": reported_site,
                "normalized_site": reported_site,
                "site_role": (
                    "matched_primary_site" if reported_site == "ovary" else "common_secondary_site"
                ),
                "sample_title": title,
                "sample_description": description,
                "source_name": source_name,
                "sra_experiment_accession": sra_matches[0],
                "pair_key_derivation": PAIR_KEY_DERIVATION,
                "metadata_anomaly": anomaly,
            }
        )

    if len({row["geo_sample_accession"] for row in rows}) != len(rows):
        raise GSE277107Error("duplicate GEO sample accessions")
    if len({row["sample_id"] for row in rows}) != len(rows):
        raise GSE277107Error("duplicate matrix sample IDs")
    by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_pair[row["pair_id"]].append(row)
    if len(by_pair) != expected_pairs:
        raise GSE277107Error(f"expected {expected_pairs} pair keys, observed {len(by_pair)}")
    for pair_id, pair_rows in by_pair.items():
        sites = [row["normalized_site"] for row in pair_rows]
        if len(pair_rows) != 2 or sorted(sites) != ["omentum", "ovary"]:
            raise GSE277107Error(f"{pair_id}: expected one ovary and one omentum sample")
    site_order = {"ovary": 0, "omentum": 1}
    rows.sort(
        key=lambda row: (
            int(row["pair_id"].removeprefix("WM")),
            site_order[row["normalized_site"]],
        )
    )
    return rows, anomalies


def _format_number(value: float) -> str:
    return format(value, ".12g")


def _write_matrix(
    path: Path, first_column: str, rows: Iterable[tuple[str, list[float]]], columns: list[str]
) -> None:
    with deterministic_gzip_text_writer(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([first_column, *columns])
        for name, values in rows:
            writer.writerow([name, *(_format_number(value) for value in values)])


def prepare_dataset(
    *,
    source_manifest_path: Path,
    source_dir: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Validate the frozen inputs and emit paired, gene-symbol CORNETO inputs."""

    source_manifest = load_source_manifest(source_manifest_path)
    sources = verify_source_files(source_manifest, source_dir)
    metadata, anomalies = build_paired_metadata(
        sources["geo_family_soft"],
        expected_samples=int(source_manifest["expected"]["rna_samples"]),
        expected_pairs=int(source_manifest["expected"]["matched_pairs"]),
    )
    metadata_by_sample = {row["sample_id"]: row for row in metadata}
    ordered_samples = [row["sample_id"] for row in metadata]

    symbol_values: dict[str, list[float]] = {}
    symbol_gene_ids: dict[str, list[tuple[str, str]]] = defaultdict(list)
    source_gene_rows = 0
    with _open_text(sources["processed_tpm"]) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise GSE277107Error("processed TPM matrix is empty") from exc
        if not header or header[0] != "GeneID":
            raise GSE277107Error("processed TPM matrix must start with GeneID")
        matrix_samples = header[1:]
        if len(matrix_samples) != len(set(matrix_samples)):
            raise GSE277107Error("processed TPM matrix contains duplicate sample columns")
        if set(matrix_samples) != set(metadata_by_sample):
            raise GSE277107Error(
                "TPM/GEO sample mismatch: "
                f"metadata_only={sorted(set(metadata_by_sample) - set(matrix_samples))}; "
                f"matrix_only={sorted(set(matrix_samples) - set(metadata_by_sample))}"
            )
        input_indices = [matrix_samples.index(sample) + 1 for sample in ordered_samples]
        seen_gene_ids: set[str] = set()
        for line_number, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise GSE277107Error(
                    f"TPM line {line_number}: {len(row)} fields, expected {len(header)}"
                )
            tokens = row[0].split("|")
            if len(tokens) != 9:
                raise GSE277107Error(f"TPM line {line_number}: unexpected GeneID encoding")
            gene_id, gene_name, biotype = tokens[0], tokens[7], tokens[8]
            if not gene_id or not gene_name or not biotype:
                raise GSE277107Error(f"TPM line {line_number}: missing gene annotation")
            if gene_id in seen_gene_ids:
                raise GSE277107Error(f"TPM line {line_number}: duplicate gene ID {gene_id}")
            seen_gene_ids.add(gene_id)
            try:
                values = [float(row[index]) for index in input_indices]
            except ValueError as exc:
                raise GSE277107Error(f"TPM line {line_number}: non-numeric value") from exc
            if any(not math.isfinite(value) or value < 0 for value in values):
                raise GSE277107Error(f"TPM line {line_number}: TPM must be finite and non-negative")
            if gene_name not in symbol_values:
                symbol_values[gene_name] = [0.0] * len(values)
            for index, value in enumerate(values):
                symbol_values[gene_name][index] += value
            symbol_gene_ids[gene_name].append((gene_id, biotype))
            source_gene_rows += 1

    expected_rows = int(source_manifest["expected"]["source_gene_rows"])
    expected_symbols = int(source_manifest["expected"]["unique_gene_symbols"])
    if source_gene_rows != expected_rows or len(symbol_values) != expected_symbols:
        raise GSE277107Error(
            "frozen expression dimensions changed: "
            f"genes={source_gene_rows}/{expected_rows}, "
            f"symbols={len(symbol_values)}/{expected_symbols}"
        )

    by_pair_site = {
        (row["pair_id"], row["normalized_site"]): ordered_samples.index(row["sample_id"])
        for row in metadata
    }
    pair_ids = sorted(
        {row["pair_id"] for row in metadata}, key=lambda item: int(item.removeprefix("WM"))
    )
    fields = [
        "study_accession",
        "geo_sample_accession",
        "sample_id",
        "pair_id",
        "aliquot_key",
        "reported_site_code",
        "reported_site",
        "normalized_site",
        "site_role",
        "sample_title",
        "sample_description",
        "source_name",
        "sra_experiment_accession",
        "pair_key_derivation",
        "metadata_anomaly",
    ]

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output_dir.parent) as temporary_name:
        stage = Path(temporary_name)
        sample_manifest = stage / "paired_sample_manifest.tsv"
        symbol_matrix = stage / "gene_symbol_tpm.tsv.gz"
        paired_delta = stage / "paired_log2_tpm_delta_omentum_minus_ovary.tsv.gz"
        gene_map = stage / "gene_id_to_symbol.tsv.gz"
        write_tsv(sample_manifest, metadata, fields)
        _write_matrix(
            symbol_matrix,
            "gene_name",
            ((gene, symbol_values[gene]) for gene in sorted(symbol_values)),
            ordered_samples,
        )
        with deterministic_gzip_text_writer(gene_map) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", "gene_name", "gene_biotype", "symbol_gene_count"])
            for gene_name in sorted(symbol_gene_ids):
                count = len(symbol_gene_ids[gene_name])
                for gene_id, biotype in sorted(symbol_gene_ids[gene_name]):
                    writer.writerow([gene_id, gene_name, biotype, count])
        _write_matrix(
            paired_delta,
            "gene_name",
            (
                (
                    gene,
                    [
                        math.log2(
                            symbol_values[gene][by_pair_site[(pair_id, "omentum")]] + 1.0
                        )
                        - math.log2(
                            symbol_values[gene][by_pair_site[(pair_id, "ovary")]] + 1.0
                        )
                        for pair_id in pair_ids
                    ],
                )
                for gene in sorted(symbol_values)
            ),
            pair_ids,
        )
        output_paths = [sample_manifest, symbol_matrix, paired_delta, gene_map]
        receipt = {
            "schema_version": "gse277107_preparation_receipt.v1",
            "status": "completed",
            "created_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
            "study_accession": "GSE277107",
            "source_manifest": {
                "path": str(source_manifest_path),
                "sha256": sha256(source_manifest_path),
            },
            "sources": {
                role: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
                for role, path in sorted(sources.items())
            },
            "validated_dimensions": {
                "rna_samples": len(metadata),
                "matched_pairs": len(pair_ids),
                "source_gene_rows": source_gene_rows,
                "unique_gene_symbols": len(symbol_values),
                "ovary_samples": sum(row["normalized_site"] == "ovary" for row in metadata),
                "omentum_samples": sum(
                    row["normalized_site"] == "omentum" for row in metadata
                ),
            },
            "derivations": {
                "pair_key": PAIR_KEY_DERIVATION,
                "gene_symbol_aggregation": "sum TPM over Ensembl rows sharing gene_name",
                "paired_delta": "log2(TPM_symbol_sum + 1) omentum minus ovary within pair",
                "clinical_patient_id_available": False,
            },
            "metadata_anomalies_preserved": anomalies,
            "claim_limits": [
                "pair_id is a metadata-derived pair key, not a clinical patient identifier",
                "bulk tissue expression includes malignant and microenvironment compartments",
                "the paired delta does not by itself establish a tumor-cell-intrinsic mechanism",
                "PRIDE PXD042150 is a linked proteomics reference but is not joined here",
            ],
            "outputs": {
                path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                for path in output_paths
            },
        }
        receipt_path = stage / "receipt.json"
        write_json(receipt_path, receipt)
        output_dir.mkdir(parents=True, exist_ok=True)
        for path in [*output_paths, receipt_path]:
            path.replace(output_dir / path.name)
    return receipt
