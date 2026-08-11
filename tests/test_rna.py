import hashlib
from pathlib import Path

from hgsoc_corneto.rna import (
    FastqSpec,
    load_rna_run_specs,
    resumable_curl_command,
    validate_fastq_file,
)
from scripts import run_salmon_quant as salmon_script

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/processed/metadata/rna_runs.tsv"


def test_all_rna_run_specs_are_paired_and_unique() -> None:
    specs = load_rna_run_specs(MANIFEST)
    assert len(specs) == 117
    assert len({spec.run_accession for spec in specs}) == 117
    assert all(len(spec.fastqs) == 2 for spec in specs)
    assert sum(fastq.bytes for spec in specs for fastq in spec.fastqs) == 477_762_645_114


def test_etab_14568_run_specs() -> None:
    specs = load_rna_run_specs(MANIFEST, study_accession="E-MTAB-14568")
    assert len(specs) == 33
    smallest = min(specs, key=lambda spec: sum(fastq.bytes for fastq in spec.fastqs))
    assert smallest.run_accession == "ERR13907062"
    assert smallest.canonical_ocm_id == "OCM341-1"
    assert smallest.fastqs[0].url.endswith("ERR13907062_1.fastq.gz")
    assert smallest.fastqs[1].url.endswith("ERR13907062_2.fastq.gz")


def test_fastq_file_validation_checks_size_and_md5(tmp_path: Path) -> None:
    target = tmp_path / "mate.fastq.gz"
    target.write_bytes(b"test payload")
    expected = FastqSpec(
        mate=1,
        url="https://example.invalid/mate.fastq.gz",
        md5="c737a42e8172ef241a45e18857b8e544",
        bytes=12,
    )
    assert validate_fastq_file(target, expected) == (True, "verified")

    wrong_size = FastqSpec(mate=1, url=expected.url, md5=expected.md5, bytes=13)
    assert validate_fastq_file(target, wrong_size) == (False, "size_mismatch:12")


def test_resumable_curl_retries_all_errors_and_stalls(tmp_path: Path) -> None:
    target = tmp_path / "run.fastq.gz.partial"
    command = resumable_curl_command(
        url="https://example.invalid/run.fastq.gz",
        target=target,
    )
    assert command[:3] == ["curl", "--location", "--fail"]
    assert command[command.index("--retry") + 1] == "12"
    assert "--retry-all-errors" in command
    assert command[command.index("--speed-limit") + 1] == "1024"
    assert command[command.index("--speed-time") + 1] == "120"
    assert command[command.index("--continue-at") + 1] == "-"
    assert command[-3:] == ["--output", str(target), "https://example.invalid/run.fastq.gz"]


def test_download_quarantines_full_checksum_mismatch(
    tmp_path: Path,
    monkeypatch,
) -> None:
    payload = b"test payload"
    target = tmp_path / "mate.fastq.gz"
    partial = target.with_name(target.name + ".partial")
    partial.write_bytes(b"fail payload")
    expected = FastqSpec(
        mate=1,
        url="https://example.invalid/mate.fastq.gz",
        md5=hashlib.md5(payload, usedforsecurity=False).hexdigest(),
        bytes=len(payload),
    )

    def fake_run(command: list[str], *, capture: bool = False) -> str:
        assert not capture
        assert command[-2] == str(partial)
        assert not partial.exists()
        partial.write_bytes(payload)
        return ""

    monkeypatch.setattr(salmon_script, "_run", fake_run)
    receipt = salmon_script._download_fastq(expected, target)

    quarantined = list(tmp_path.glob("mate.fastq.gz.partial.invalid-*"))
    assert target.read_bytes() == payload
    assert len(quarantined) == 1
    assert quarantined[0].read_bytes() == b"fail payload"
    assert receipt["verification"] == "verified"
    assert receipt["quarantined_partials"][0]["reason"].startswith("md5_mismatch:")
