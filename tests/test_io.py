from pathlib import Path

from hgsoc_corneto.io import deterministic_gzip_text_writer, sha256


def test_deterministic_gzip_writer(tmp_path: Path) -> None:
    paths = (tmp_path / "first.tsv.gz", tmp_path / "second.tsv.gz")
    for path in paths:
        with deterministic_gzip_text_writer(path) as handle:
            handle.write("a\tb\n1\t2\n")
    assert sha256(paths[0]) == sha256(paths[1])
