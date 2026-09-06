import json

from audit_fixed_indicator_feasibility import restricted_lp
from audit_metabolic_scientific_validity import bound_violations, inspect_solution


def example_model():
    return {
        "uptake": {"lower": 0, "upper": 10, "stoichiometry": {"A": 1}},
        "biomass_human": {"lower": 0, "upper": 100, "stoichiometry": {"A": -1}},
    }


def test_nonselected_uptake_cannot_sustain_required_growth():
    model = example_model()
    assert (
        restricted_lp(model, {"biomass_human": (9, 100)}, {"biomass_human"})["status"]
        == "infeasible"
    )
    result = restricted_lp(model, {"biomass_human": (9, 100)}, set(model))
    assert result["success"]
    assert abs(result["maximum_biomass"] - 10) < 1e-8


def test_saved_growth_is_a_separate_feasibility_test():
    model = example_model()
    result = restricted_lp(model, {"biomass_human": (9, 100)}, set(model), required_biomass=11)
    assert result["status"] == "infeasible"


def test_sparse_diagnostics_check_numerics_but_do_not_imply_specificity():
    solution = {
        "nonzero_fluxes": [["uptake", 1], ["biomass_human", 1]],
        "active_by_indicator": ["biomass_human"],
    }
    result = inspect_solution(
        solution, example_model(), {"biomass_human": [0.9, 100]}, {"biomass_human": -1}
    )
    assert result["flux_without_selected_indicator_count"] == 1
    assert result["max_sparse_summary_mass_balance_residual"] == 0
    assert result["expression_cap_count"] == 0
    assert bound_violations(dict(solution["nonzero_fluxes"]), {"uptake": [0, 0.5]}) == ["uptake"]
    json.dumps(result, allow_nan=False)
