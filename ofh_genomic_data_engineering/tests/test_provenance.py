from pathlib import Path

from scripts.build_provenance import sha256


def test_sha256_is_deterministic(tmp_path: Path) -> None:
    p = tmp_path / "x"
    p.write_bytes(b"abc")
    assert sha256(p) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
