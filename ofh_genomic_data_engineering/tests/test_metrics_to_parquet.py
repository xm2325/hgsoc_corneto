from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from scripts.metrics_to_parquet import build_outputs, read_table, semantic_table_sha256


def write(path: Path, text: str) -> None:
    path.write_text(text)


def test_read_table_strips_hash(tmp_path: Path) -> None:
    p = tmp_path / "x.tsv"
    write(p, "#CHROM\tPOS\tID\n22\t1\trs1\n")
    df = read_table(p)
    assert list(df.columns) == ["CHROM", "POS", "ID"]
    assert df.iloc[0]["ID"] == "rs1"


def test_read_table_skips_vcf_style_metadata(tmp_path: Path) -> None:
    p = tmp_path / "x.pvar"
    write(p, "##fileformat=VCFv4.2\n##contig=<ID=22>\n#CHROM\tPOS\tID\tREF\tALT\n22\t16050075\trs587697622\tA\tG\n")
    df = read_table(p)
    assert list(df.columns) == ["CHROM", "POS", "ID", "REF", "ALT"]
    assert len(df) == 1
    assert df.iloc[0]["POS"] == "16050075"


def test_read_table_rejects_missing_header(tmp_path: Path) -> None:
    p = tmp_path / "bad.txt"
    write(p, "plain text without a PLINK header\n")
    with pytest.raises(ValueError, match="no PLINK table header"):
        read_table(p)


def test_build_outputs_writes_typed_tables_and_schema(tmp_path: Path) -> None:
    pytest.importorskip("pyarrow")
    files = {}
    contents = {
        "pvar": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n22\t2\trs2\tC\tT\n22\t1\trs1\tA\tG\n",
        "psam": "#FID\tIID\nS1\tS1\nS2\tS2\n",
        "afreq": "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n22\trs1\tA\tG\t0.25\t4\n",
        "vmiss": "#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS\n22\trs1\t0\t2\t0\n",
        "smiss": "#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS\nS1\tS1\t0\t2\t0\n",
        "hardy": "#CHROM\tID\tA1\tAX\tHOM_A1_CT\tHET_A1_CT\tTWO_AX_CT\tO_HET\tE_HET\tP\n22\trs1\tG\tA\t0\t1\t1\t0.5\t0.5\t1\n",
        "eigenvec": "#FID\tIID\tPC1\tPC2\nS1\tS1\t0.1\t0.2\n",
    }
    for key, text in contents.items():
        p = tmp_path / key
        write(p, text)
        files[key] = str(p)
    schema_manifest = tmp_path / "schema_manifest.json"
    args = Namespace(**files, outdir=str(tmp_path / "pq"), summary=str(tmp_path / "summary.json"), schema_manifest=str(schema_manifest))
    summary = build_outputs(args)
    assert summary["sample_count"] == 2
    assert summary["variant_count"] == 2
    assert summary["storage"] == {"format": "parquet", "compression": "zstd"}
    assert len(list((tmp_path / "pq").glob("*.parquet"))) == 7
    assert len(summary["semantic_hashes"]) == 7
    assert all(len(value) == 64 for value in summary["semantic_hashes"].values())
    variants = pd.read_parquet(tmp_path / "pq" / "variants.parquet")
    assert variants["ID"].tolist() == ["rs1", "rs2"]
    assert str(variants["POS"].dtype) == "int64"
    assert summary["semantic_hashes"]["variants"] == semantic_table_sha256(variants)
    afreq = pd.read_parquet(tmp_path / "pq" / "allele_frequencies.parquet")
    assert str(afreq["ALT_FREQS"].dtype) == "float64"
    assert str(afreq["OBS_CT"].dtype) == "int64"
    import json
    manifest = json.loads(schema_manifest.read_text())
    variant_types = {c["name"]: c["type"] for c in manifest["tables"]["variants"]["columns"]}
    assert variant_types["POS"] == "int64"


def test_semantic_hash_changes_with_logical_content() -> None:
    original = pd.DataFrame({"IID": ["S1", "S2"], "PC1": [0.1, -0.1]})
    changed = original.copy()
    changed.loc[0, "PC1"] = 0.2
    assert semantic_table_sha256(original) != semantic_table_sha256(changed)
