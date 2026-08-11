#!/usr/bin/env python3
"""Fetch the small public metadata needed to rebuild the OCM manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STUDIES = {
    "E-MTAB-7223": "PRJEB28709",
    "E-MTAB-10801": "PRJEB46736",
    "E-MTAB-11000": "PRJEB47842",
    "E-MTAB-14568": "PRJEB81794",
}
ENA_FIELDS = ",".join(
    [
        "study_accession",
        "secondary_study_accession",
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "secondary_sample_accession",
        "sample_alias",
        "sample_title",
        "library_name",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
        "first_public",
    ]
)


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as handle:
        temporary = Path(handle.name)
        request = urllib.request.Request(url, headers={"User-Agent": "hgsoc-corneto/0.1"})
        with urllib.request.urlopen(request, timeout=120) as response:
            shutil.copyfileobj(response, handle)
    temporary.replace(target)


def md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def fetch_study(study: str, project: str, output: Path) -> None:
    api = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{study}"
    info_path = output / "biostudies" / f"{study}.info.json"
    download(f"{api}/info", info_path)
    download(api, output / "biostudies" / f"{study}.study.json")
    info = json.loads(info_path.read_text(encoding="utf-8"))
    root = info["httpLink"]
    for kind in ("idf", "sdrf"):
        name = f"{study}.{kind}.txt"
        download(f"{root}/Files/{name}", output / "biostudies" / name)

    query = urllib.parse.urlencode(
        {
            "accession": project,
            "result": "read_run",
            "fields": ENA_FIELDS,
            "format": "tsv",
        }
    )
    download(
        f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}",
        output / "ena" / f"{project}.read_run.tsv",
    )


def fetch_tighe(output: Path, supplement_dir: Path) -> None:
    pmcid = "PMC12208324"
    download(
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
        output / f"{pmcid}.fullText.xml",
    )
    zip_path = supplement_dir / f"{pmcid}.supplementaryFiles.zip"
    download(
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles",
        zip_path,
    )
    with zipfile.ZipFile(zip_path) as archive:
        for member, destination, expected_md5 in [
            ("mmc1.pdf", supplement_dir / "mmc1.pdf", "dd534b92b1b8c64bcba970127ee16767"),
            ("mmc2.xlsx", output / "mmc2.xlsx", "fc47d624e4d7ef1ae87394357189b049"),
        ]:
            destination.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, destination.open("wb") as target:
                shutil.copyfileobj(source, target)
            if md5(destination) != expected_md5:
                raise ValueError(f"MD5 mismatch for {member}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "data/raw/metadata")
    parser.add_argument("--supplement-dir", type=Path, default=ROOT / "tmp/pdfs")
    parser.add_argument("--skip-supplement", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for study, project in STUDIES.items():
        fetch_study(study, project, args.output)
    if not args.skip_supplement:
        fetch_tighe(args.output, args.supplement_dir)


if __name__ == "__main__":
    main()
