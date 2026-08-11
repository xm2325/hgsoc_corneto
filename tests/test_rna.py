from pathlib import Path

from hgsoc_corneto.rna import FastqSpec, load_rna_run_specs, validate_fastq_file

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
