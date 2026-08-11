from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "build_pooled_primary_expression.py"
SPEC = importlib.util.spec_from_file_location("build_pooled_primary_expression", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_matrix(path: Path, samples: list[str], offset: int) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_name", *samples])
        writer.writerow(["ENSG1", "A", *(offset + index for index in range(len(samples)))])
        writer.writerow(["ENSG2", "B", *(offset + 10 + index for index in range(len(samples)))])


def test_build_pooled_matrix_selects_only_primary_runs(tmp_path: Path) -> None:
    fields = [
        "study_accession", "run_accession", "canonical_ocm_id", "patient_id",
        "sample_class", "histotype_group", "is_representative_rna_library",
        "primary_cohort_eligible", "chemo_naive_at_biopsy", "biopsy_type",
    ]
    rows = [
        ["S1", "R1", "OCM1", "P1", "tumour", "HGSOC", "true", "true", "true", "Ascites"],
        ["S1", "R2", "OCM2", "P2", "stroma", "HGSOC", "true", "false", "true", "Ascites"],
        ["S2", "R3", "OCM3", "P3", "tumour", "HGSOC", "true", "true", "false", "Solid"],
    ]
    manifest = tmp_path / "manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    first = tmp_path / "s1.tsv.gz"
    second = tmp_path / "s2.tsv.gz"
    _write_matrix(first, ["R2", "R1"], 0)
    _write_matrix(second, ["R3"], 20)
    output = tmp_path / "out"
    receipt = MODULE.build_pooled_matrix(
        manifest_path=manifest,
        matrices=[("S1", first), ("S2", second)],
        output_dir=output,
        value_name="tpm",
        expected_samples=2,
        expected_patients=2,
    )
    assert receipt["primary_run_count"] == 2
    assert receipt["study_counts"] == {"S1": 1, "S2": 1}
    with gzip.open(output / "gene_tpm.tsv.gz", "rt", encoding="utf-8") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    assert table[0] == ["gene_id", "gene_name", "R1", "R3"]
    assert table[1] == ["ENSG1", "A", "1", "20"]
