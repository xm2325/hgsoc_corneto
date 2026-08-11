from pathlib import Path

import pytest

from hgsoc_corneto.rna import (
    aggregate_salmon_quant,
    iter_gene_matrix_rows,
    load_gencode_gene_map,
    parse_gtf_attributes,
)


def test_parse_gtf_attributes() -> None:
    assert parse_gtf_attributes('gene_id "ENSG1.1"; gene_name "A"; tag "basic";') == {
        "gene_id": "ENSG1.1",
        "gene_name": "A",
        "tag": "basic",
    }


def test_salmon_transcripts_aggregate_to_genes(tmp_path: Path) -> None:
    gtf = tmp_path / "toy.gtf"
    gtf.write_text(
        "##gtf-version 2.2\n"
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG1.1"; '
        'transcript_id "ENST1.1"; gene_type "protein_coding"; gene_name "A";\n'
        'chr1\ttest\ttranscript\t20\t90\t.\t+\t.\tgene_id "ENSG1.1"; '
        'transcript_id "ENST2.1"; gene_type "protein_coding"; gene_name "A";\n'
        'chr2\ttest\ttranscript\t5\t50\t.\t-\t.\tgene_id "ENSG2.2"; '
        'transcript_id "ENST3.2"; gene_type "lncRNA"; gene_name "B";\n',
        encoding="utf-8",
    )
    quant = tmp_path / "quant.sf"
    quant.write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "ENST1.1|ENSG1.1\t100\t80\t100000\t5.5\n"
        "ENST2.1\t80\t60\t300000\t7.25\n"
        "ENST3.2\t50\t30\t600000\t9\n",
        encoding="utf-8",
    )
    genes, transcript_map = load_gencode_gene_map(gtf)
    assert [gene.gene_id for gene in genes] == ["ENSG1.1", "ENSG2.2"]
    assert genes[0].transcript_count == 2
    assert genes[0].start == 1
    assert genes[0].end == 100

    sample = aggregate_salmon_quant(
        quant,
        run_accession="ERR1",
        transcript_to_gene_index=transcript_map,
        gene_count=len(genes),
    )
    assert sample.counts == pytest.approx((12.75, 9.0))
    assert sample.tpm == pytest.approx((400000, 600000))
    assert sample.transcript_rows == sample.mapped_transcript_rows == 3
    assert sample.unmapped_transcript_ids == ()
    rows = list(iter_gene_matrix_rows(genes, (sample,), value="counts"))
    assert rows[0] == ("ENSG1.1", "A", (12.75,))


def test_unmapped_transcript_is_auditable(tmp_path: Path) -> None:
    quant = tmp_path / "quant.sf"
    quant.write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "UNKNOWN.1\t100\t80\t1\t2\n",
        encoding="utf-8",
    )
    sample = aggregate_salmon_quant(
        quant,
        run_accession="ERR2",
        transcript_to_gene_index={},
        gene_count=1,
    )
    assert sample.unmapped_transcript_ids == ("UNKNOWN.1",)
