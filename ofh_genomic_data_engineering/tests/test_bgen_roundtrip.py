from pathlib import Path

from scripts.validate_bgen_roundtrip import validate_bgen_roundtrip


def _write_fixture(tmp_path: Path) -> dict[str, Path]:
    source_pvar = tmp_path / "source.pvar"
    source_psam = tmp_path / "source.psam"
    source_afreq = tmp_path / "source.afreq"
    roundtrip_pvar = tmp_path / "roundtrip.pvar"
    roundtrip_psam = tmp_path / "roundtrip.psam"
    roundtrip_afreq = tmp_path / "roundtrip.afreq"

    pvar = "#CHROM POS ID REF ALT\n22 101 rs1 A G\n22 202 rs2 C T\n"
    psam = "#FID IID\nS1 S1\nS2 S2\n"
    afreq = "#CHROM ID REF ALT ALT_FREQS OBS_CT\n22 rs1 A G 0.25 4\n22 rs2 C T 0.5 4\n"

    source_pvar.write_text(pvar)
    source_psam.write_text(psam)
    source_afreq.write_text(afreq)
    roundtrip_pvar.write_text(pvar)
    roundtrip_psam.write_text(psam)
    roundtrip_afreq.write_text(afreq)

    return {
        "source_pvar": source_pvar,
        "source_psam": source_psam,
        "source_afreq": source_afreq,
        "roundtrip_pvar": roundtrip_pvar,
        "roundtrip_psam": roundtrip_psam,
        "roundtrip_afreq": roundtrip_afreq,
        "output_path": tmp_path / "bgen_validation.json",
    }


def test_bgen_roundtrip_passes_for_matching_semantics(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    payload = validate_bgen_roundtrip(**paths, frequency_tolerance=1e-4)
    assert payload["status"] == "PASS"
    assert payload["allele_convention"] == "ref-first"
    assert payload["probability_bits"] == 16
    assert payload["sample_count"] == 2
    assert payload["variant_count"] == 2
    assert payload["max_abs_alt_frequency_diff"] == 0.0
    assert all(check["status"] == "PASS" for check in payload["checks"])


def test_bgen_roundtrip_rejects_sample_order_change(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    paths["roundtrip_psam"].write_text("#FID IID\nS2 S2\nS1 S1\n")
    payload = validate_bgen_roundtrip(**paths, frequency_tolerance=1e-4)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "sample_order_identity" in failed
    assert "sample_identity_hash" in failed


def test_bgen_roundtrip_rejects_ref_alt_flip(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    paths["roundtrip_pvar"].write_text("#CHROM POS ID REF ALT\n22 101 rs1 G A\n22 202 rs2 C T\n")
    paths["roundtrip_afreq"].write_text(
        "#CHROM ID REF ALT ALT_FREQS OBS_CT\n22 rs1 G A 0.75 4\n22 rs2 C T 0.5 4\n"
    )
    payload = validate_bgen_roundtrip(**paths, frequency_tolerance=1e-4)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "variant_order_allele_identity" in failed
    assert "frequency_variant_order_identity" in failed


def test_bgen_roundtrip_rejects_frequency_drift(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    paths["roundtrip_afreq"].write_text(
        "#CHROM ID REF ALT ALT_FREQS OBS_CT\n22 rs1 A G 0.251 4\n22 rs2 C T 0.5 4\n"
    )
    payload = validate_bgen_roundtrip(**paths, frequency_tolerance=1e-4)
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "alt_frequency_tolerance" in failed
