from pathlib import Path

import numpy as np

from scripts.run_nmf_primary import load_primary_cohorts


def _write_inputs(root: Path, study: str, run_prefix: str) -> tuple[str, Path, Path]:
    matrix = root / f"{study}.tsv"
    qc = root / f"{study}.qc.tsv"
    runs = [f"{run_prefix}1", f"{run_prefix}2", f"{run_prefix}3"]
    matrix.write_text(
        "gene_id\tgene_name\t" + "\t".join(runs) + "\n"
        "G1\tA\t1\t2\t3\n"
        "G2\tB\t3\t2\t1\n",
        encoding="utf-8",
    )
    qc.write_text(
        "study_accession\trun_accession\tcanonical_ocm_id\tpatient_id\tsample_class\thistotype_group\tprimary_cohort_eligible\n"
        f"{study}\t{runs[0]}\tOCM1\tP1\ttumour\tHGSOC\ttrue\n"
        f"{study}\t{runs[1]}\tOCM2\tP2\ttumour\tHGSOC\ttrue\n"
        f"{study}\t{runs[2]}\tNA\tNA\tcell_line_control\tNA\tfalse\n",
        encoding="utf-8",
    )
    return study, matrix, qc


def test_load_primary_cohorts_filters_and_preserves_study(tmp_path: Path) -> None:
    first = _write_inputs(tmp_path, "E-MTAB-A", "ERR1")
    second = _write_inputs(tmp_path, "E-MTAB-B", "ERR2")
    genes, names, samples, values, metadata, inputs = load_primary_cohorts([first, second])
    assert genes == ("G1", "G2")
    assert names == ("A", "B")
    assert samples == ("ERR11", "ERR12", "ERR21", "ERR22")
    assert values.shape == (2, 4)
    assert [row["study_accession"] for row in metadata] == ["E-MTAB-A"] * 2 + ["E-MTAB-B"] * 2
    assert [row["primary_samples"] for row in inputs] == [2, 2]
    assert np.all(values >= 0)
