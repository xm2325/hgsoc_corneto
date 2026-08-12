import csv
import gzip
import hashlib
import json
from pathlib import Path

from scripts.prepare_gse208216_corneto_input import prepare


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_prepare_symbol_cpm_and_model_manifest(tmp_path: Path) -> None:
    counts = tmp_path / "counts.tsv.gz"
    with gzip.open(counts, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(
            [
                ["", "PDO1", "FT1"],
                ["ENSG1.1", 1, 2],
                ["ENSG2.2", 3, 2],
                ["ENSG3.1", 1, 1],
            ]
        )
    gene_map = tmp_path / "gene_metadata.tsv"
    gene_map.write_text(
        "gene_id\tgene_name\nENSG1.9\tA\nENSG2.2\tA\nENSG3.1\tB\n",
        encoding="utf-8",
    )
    contract = tmp_path / "contract.json"
    contract.write_text(
        json.dumps(
            {
                "count_matrix": {"sha256": _sha(counts), "gene_rows": 3},
                "samples": [
                    {
                        "sample_id": "PDO1",
                        "sample_class": "hgsoc_organoid",
                        "geo_accession": "GSM1",
                    },
                    {
                        "sample_id": "FT1",
                        "sample_class": "fallopian_tube_organoid",
                        "geo_accession": "GSM2",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "out"
    result = prepare(counts, gene_map, contract, output)
    assert result["status"] == "completed"
    assert result["dimensions"]["unique_gene_symbols"] == 2
    with gzip.open(
        output / "gene_symbol_cpm.tsv.gz", "rt", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0] == ["gene_name", "PDO1", "FT1"]
    assert float(rows[1][1]) == 800_000
    assert float(rows[1][2]) == 800_000
    manifest = (output / "corneto_manifest.tsv").read_text(encoding="utf-8")
    assert "hgsoc_organoid" in manifest and "fallopian_tube_organoid" in manifest
