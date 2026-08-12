import csv
import gzip
import json
from pathlib import Path

import pytest

from scripts.normalize_external_pseudobulk_cpm import normalize


def _write(path: Path, rows: list[list[object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def test_normalize_raw_counts_to_dataset_internal_cpm(tmp_path: Path) -> None:
    source = tmp_path / "counts.tsv.gz"
    _write(source, [["gene", "S1", "S2"], ["A", 1, 3], ["B", 1, 1]])
    output = tmp_path / "cpm.tsv.gz"
    receipt = tmp_path / "receipt.json"
    result = normalize(source, output, receipt, "TEST", 2, 2)
    assert result["status"] == "completed"
    with gzip.open(output, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0] == ["gene_name", "S1", "S2"]
    assert float(rows[1][1]) == 500_000
    assert float(rows[1][2]) == 750_000
    assert json.loads(receipt.read_text())["dimensions"] == {"genes": 2, "samples": 2}


def test_zero_library_fails_closed(tmp_path: Path) -> None:
    source = tmp_path / "counts.tsv.gz"
    _write(source, [["gene", "S1", "S2"], ["A", 1, 0]])
    with pytest.raises(ValueError, match="zero-size library"):
        normalize(source, tmp_path / "out.gz", tmp_path / "r.json", "TEST", 2, 1)
