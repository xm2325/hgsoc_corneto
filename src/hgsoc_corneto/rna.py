"""Validated run specifications for restartable Salmon quantification."""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from pathlib import Path
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
