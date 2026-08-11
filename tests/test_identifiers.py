from hgsoc_corneto.identifiers import (
    apply_alias,
    canonicalize_tighe_ocm_id,
    parse_source_name,
)


def test_parse_legacy_punctuation_and_passage():
    parsed = parse_source_name("OCM.46-3T-P14")
    assert parsed.source_biospecimen_id == "OCM46-3"
    assert parsed.patient_id == "OCM46"
    assert parsed.material == "tumour"
    assert parsed.passage == "P14"


def test_parse_fractionated_ocm():
    parsed = parse_source_name("OCM.64-3-T-EpcamNeg")
    assert parsed.source_biospecimen_id == "OCM64-3-EPCAMNEG"
    assert parsed.material == "tumour"


def test_parse_alpha_spatial_model():
    parsed = parse_source_name("OCM361a-T")
    assert parsed.source_biospecimen_id == "OCM361a"
    assert parsed.patient_id == "OCM361"


def test_controls_are_not_ocms():
    assert parse_source_name("Kuramochi").is_control
    assert parse_source_name("FNE").source_biospecimen_id is None


def test_explicit_alias_is_applied_without_changing_patient():
    assert apply_alias("OCM231", {"OCM231": "OCM231-1"}) == "OCM231-1"


def test_tighe_table_s2_notation():
    assert canonicalize_tighe_ocm_id("64-3-") == "OCM64-3-Ep-"
    assert canonicalize_tighe_ocm_id("376Ta") == "OCM376a"
