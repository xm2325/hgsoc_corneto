import json

import pytest
from validate_revised_metabolic_model import validate_support


def test_fixed_support_and_json():
    reactions = {
        "uptake": {"lower": 0, "upper": 10, "stoichiometry": {"a": 1}},
        "convert": {"lower": 0, "upper": 100, "stoichiometry": {"a": -1, "b": 1}},
        "biomass_human": {"lower": 0, "upper": 100, "stoichiometry": {"b": -1}},
        "unused": {"lower": -100, "upper": 100, "stoichiometry": {"c": 1}},
    }
    r = validate_support(reactions, {"convert": (0, 4)}, fraction=0.9)
    assert r["status"] == "validated_support"
    assert r["pfba_biomass"] == pytest.approx(3.6)
    assert "unused" not in r["selected_reactions"]
    assert r["fixed_support_check"]["maximum_biomass"] == pytest.approx(4)
    json.dumps(r, allow_nan=False)
    r = validate_support(reactions, {"uptake": (0, 0)})
    assert r["status"] == "no_positive_growth"
