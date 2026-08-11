"""Validated run specifications and Salmon-to-gene aggregation utilities."""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO
from urllib.parse import urlparse

from hgsoc_corneto.io import read_tsv


@dataclass(frozen=True)
class FastqSpec:
    mate: int
    url: str
    md5: str
    bytes: int


@dataclass(frozen=True)
class RnaRunSpec:
    study_accession: str
    run_accession: str
    canonical_ocm_id: str | None
    patient_id: str | None
    sample_class: str
    library_layout: str
    fastqs: tuple[FastqSpec, FastqSpec]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class GeneRecord:
    gene_id: str
    gene_name: str
    gene_type: str
    chromosome: str
    start: int
    end: int
    strand: str
    transcript_count: int


@dataclass(frozen=True)
class SalmonGeneSample:
    run_accession: str
    counts: tuple[float, ...]
    tpm: tuple[float, ...]
    transcript_rows: int
    mapped_transcript_rows: int
    unmapped_transcript_ids: tuple[str, ...]

    @property
    def estimated_count_sum(self) -> float:
        return sum(self.counts)

    @property
    def tpm_sum(self) -> float:
        return sum(self.tpm)


def _optional(value: str) -> str | None:
    return None if value in {"", "NA"} else value


def _https_url(value: str) -> str:
    value = value.strip()
    if value.startswith("ftp://"):
        value = "https://" + value.removeprefix("ftp://")
    elif "://" not in value:
        value = "https://" + value
    parsed = urlparse(value)
    if parsed.scheme != "https" or not parsed.netloc or not parsed.path:
        raise ValueError(f"Invalid ENA FASTQ URL: {value!r}")
    return value


def _fastq_specs(row: dict[str, str]) -> tuple[FastqSpec, FastqSpec]:
    urls = row["fastq_ftp"].split(";")
    md5s = row["fastq_md5"].split(";")
    sizes = row["fastq_bytes"].split(";")
    if row["library_layout"] != "PAIRED":
        raise ValueError(f"Only paired libraries are supported: {row['run_accession']}")
    if not (len(urls) == len(md5s) == len(sizes) == 2):
        raise ValueError(f"Expected two FASTQs for {row['run_accession']}")
    specs = tuple(
        FastqSpec(
            mate=mate,
            url=_https_url(url),
            md5=checksum.lower(),
            bytes=int(size),
        )
        for mate, (url, checksum, size) in enumerate(
            zip(urls, md5s, sizes, strict=True), start=1
        )
    )
    for spec in specs:
        if len(spec.md5) != 32 or any(char not in "0123456789abcdef" for char in spec.md5):
            raise ValueError(f"Invalid MD5 for {row['run_accession']} mate {spec.mate}")
        if spec.bytes <= 0:
            raise ValueError(f"Invalid FASTQ size for {row['run_accession']} mate {spec.mate}")
    return specs  # type: ignore[return-value]


def load_rna_run_specs(
    manifest: str | Path,
    *,
    study_accession: str | None = None,
) -> tuple[RnaRunSpec, ...]:
    rows = read_tsv(manifest)
    if study_accession is not None:
        rows = [row for row in rows if row["study_accession"] == study_accession]
    specs = tuple(
        RnaRunSpec(
            study_accession=row["study_accession"],
            run_accession=row["run_accession"],
            canonical_ocm_id=_optional(row["canonical_ocm_id"]),
            patient_id=_optional(row["patient_id"]),
            sample_class=row["sample_class"],
            library_layout=row["library_layout"],
            fastqs=_fastq_specs(row),
        )
        for row in rows
    )
    run_ids = [spec.run_accession for spec in specs]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("Run accessions must be unique after filtering")
    if study_accession is not None and not specs:
        raise ValueError(f"No runs found for study {study_accession}")
    return specs


def file_md5(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def validate_fastq_file(path: str | Path, spec: FastqSpec) -> tuple[bool, str]:
    target = Path(path)
    if not target.is_file():
        return False, "missing"
    actual_size = target.stat().st_size
    if actual_size != spec.bytes:
        return False, f"size_mismatch:{actual_size}"
    actual_md5 = file_md5(target)
    if actual_md5 != spec.md5:
        return False, f"md5_mismatch:{actual_md5}"
    return True, "verified"


def parse_gtf_attributes(value: str) -> dict[str, str]:
    """Parse the semicolon-delimited attribute field from a GENCODE GTF row."""

    attributes: dict[str, str] = {}
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, raw_value = item.partition(" ")
        if not separator:
            raise ValueError(f"Malformed GTF attribute: {item!r}")
        parsed_value = raw_value.strip()
        if len(parsed_value) >= 2 and parsed_value[0] == parsed_value[-1] == '"':
            parsed_value = parsed_value[1:-1]
        if key in attributes and attributes[key] != parsed_value:
            raise ValueError(f"Conflicting GTF attribute {key!r}")
        attributes[key] = parsed_value
    return attributes


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def load_gencode_gene_map(
    gtf_path: str | Path,
) -> tuple[tuple[GeneRecord, ...], dict[str, int]]:
    """Return ordered GENCODE genes and a versioned transcript-to-gene index."""

    target = Path(gtf_path)
    transcript_rows: list[tuple[str, str, str, str, str, int, int, str]] = []
    gene_order: list[str] = []
    seen_genes: set[str] = set()
    with _open_text(target) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected 9 GTF columns at {target}:{line_number}")
            chromosome, _source, feature, start, end, _score, strand, _frame, raw_attrs = fields
            if feature != "transcript":
                continue
            attrs = parse_gtf_attributes(raw_attrs)
            try:
                transcript_id = attrs["transcript_id"]
                gene_id = attrs["gene_id"]
            except KeyError as error:
                raise ValueError(
                    f"Missing {error.args[0]} at {target}:{line_number}"
                ) from error
            gene_name = attrs.get("gene_name", gene_id)
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", "NA"))
            transcript_rows.append(
                (
                    transcript_id,
                    gene_id,
                    gene_name,
                    gene_type,
                    chromosome,
                    int(start),
                    int(end),
                    strand,
                )
            )
            if gene_id not in seen_genes:
                gene_order.append(gene_id)
                seen_genes.add(gene_id)
    if not transcript_rows:
        raise ValueError(f"No transcript features found in {target}")

    gene_values: dict[str, dict[str, object]] = {}
    transcript_to_gene_id: dict[str, str] = {}
    for (
        transcript_id,
        gene_id,
        gene_name,
        gene_type,
        chromosome,
        start,
        end,
        strand,
    ) in transcript_rows:
        previous_gene = transcript_to_gene_id.setdefault(transcript_id, gene_id)
        if previous_gene != gene_id:
            raise ValueError(f"Transcript maps to multiple genes: {transcript_id}")
        if gene_id not in gene_values:
            gene_values[gene_id] = {
                "gene_name": gene_name,
                "gene_type": gene_type,
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "strand": strand,
                "transcript_count": 0,
            }
        values = gene_values[gene_id]
        if values["gene_name"] != gene_name or values["gene_type"] != gene_type:
            raise ValueError(f"Conflicting metadata for gene {gene_id}")
        if values["chromosome"] != chromosome or values["strand"] != strand:
            raise ValueError(f"Conflicting coordinates for gene {gene_id}")
        values["start"] = min(int(values["start"]), start)
        values["end"] = max(int(values["end"]), end)
        values["transcript_count"] = int(values["transcript_count"]) + 1

    genes = tuple(
        GeneRecord(
            gene_id=gene_id,
            gene_name=str(gene_values[gene_id]["gene_name"]),
            gene_type=str(gene_values[gene_id]["gene_type"]),
            chromosome=str(gene_values[gene_id]["chromosome"]),
            start=int(gene_values[gene_id]["start"]),
            end=int(gene_values[gene_id]["end"]),
            strand=str(gene_values[gene_id]["strand"]),
            transcript_count=int(gene_values[gene_id]["transcript_count"]),
        )
        for gene_id in gene_order
    )
    gene_index = {gene.gene_id: index for index, gene in enumerate(genes)}
    return genes, {
        transcript_id: gene_index[gene_id]
        for transcript_id, gene_id in transcript_to_gene_id.items()
    }


def _salmon_transcript_id(value: str) -> str:
    """Normalize GENCODE FASTA headers whether Salmon kept or split pipe fields."""

    return value.split("|", maxsplit=1)[0]


def aggregate_salmon_quant(
    quant_path: str | Path,
    *,
    run_accession: str,
    transcript_to_gene_index: dict[str, int],
    gene_count: int,
) -> SalmonGeneSample:
    """Sum transcript-level Salmon NumReads and TPM values to GENCODE genes."""

    counts = [0.0] * gene_count
    tpm = [0.0] * gene_count
    transcript_rows = 0
    mapped_rows = 0
    unmapped: list[str] = []
    seen: set[str] = set()
    with Path(quant_path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Name", "TPM", "NumReads"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid Salmon quant.sf header: {quant_path}")
        for row in reader:
            transcript_rows += 1
            transcript_id = _salmon_transcript_id(row["Name"])
            if transcript_id in seen:
                raise ValueError(f"Duplicate Salmon target: {transcript_id}")
            seen.add(transcript_id)
            gene_index = transcript_to_gene_index.get(transcript_id)
            if gene_index is None:
                unmapped.append(transcript_id)
                continue
            count_value = float(row["NumReads"])
            tpm_value = float(row["TPM"])
            if not math.isfinite(count_value) or count_value < 0:
                raise ValueError(f"Invalid NumReads for {transcript_id}: {count_value}")
            if not math.isfinite(tpm_value) or tpm_value < 0:
                raise ValueError(f"Invalid TPM for {transcript_id}: {tpm_value}")
            counts[gene_index] += count_value
            tpm[gene_index] += tpm_value
            mapped_rows += 1
    if transcript_rows == 0:
        raise ValueError(f"No transcript rows in {quant_path}")
    return SalmonGeneSample(
        run_accession=run_accession,
        counts=tuple(counts),
        tpm=tuple(tpm),
        transcript_rows=transcript_rows,
        mapped_transcript_rows=mapped_rows,
        unmapped_transcript_ids=tuple(sorted(unmapped)),
    )


def iter_gene_matrix_rows(
    genes: tuple[GeneRecord, ...],
    samples: tuple[SalmonGeneSample, ...],
    *,
    value: str,
) -> Iterator[tuple[str, str, tuple[float, ...]]]:
    """Yield one gene row for counts, TPM, or log1p(TPM)."""

    if value not in {"counts", "tpm", "log1p_tpm"}:
        raise ValueError(f"Unknown matrix value: {value}")
    for index, gene in enumerate(genes):
        if value == "counts":
            values = tuple(sample.counts[index] for sample in samples)
        elif value == "tpm":
            values = tuple(sample.tpm[index] for sample in samples)
        else:
            values = tuple(math.log1p(sample.tpm[index]) for sample in samples)
        yield gene.gene_id, gene.gene_name, values
