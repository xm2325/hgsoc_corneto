import json
from types import SimpleNamespace

import numpy as np
import pytest
from interpret_existing_biology import paired_stroma
from research_validation_common import bh, heldout_split, one_per_patient, select_mad, signflip
from run_metabolic_information_pilot import (
    ContinuousModel,
    block_uptake,
    cap_policy,
    evaluate,
    normalize_gene_id,
    normalize_unique_gene_ids,
)
from run_patient_grouped_nmf_validation import fit_best, normalized_labels


def test_gene_version_normalization_preserves_identity():
    assert normalize_gene_id("ENSG000001.12") == "ENSG000001"
    assert normalize_gene_id("ENSG000001.12_PAR_Y") == "ENSG000001_PAR_Y"
    assert normalize_unique_gene_ids(["ENSG000001.12", "ENSG000001.12_PAR_Y"]) == [
        "ENSG000001", "ENSG000001_PAR_Y"
    ]
    for identifier in ["GENE.1", "ENSG000001.bad", "ENSG000001_PAR_Y", "ENSG000001.2_other"]:
        assert normalize_gene_id(identifier) == identifier


def test_gene_normalization_rejects_real_collisions():
    with pytest.raises(ValueError, match="ENSG000001"):
        normalize_unique_gene_ids(["ENSG000001.1", "ENSG000001.2"])


def test_lp_receipt_serializes_numpy_bound_comparison(monkeypatch):
    monkeypatch.setattr("run_metabolic_information_pilot.linprog", lambda *a, **kw:
        SimpleNamespace(status=0, message="test", success=True,
                        x=np.array([-1e-10]), fun=1e-10))
    model = ContinuousModel({"biomass_human": {
        "lower": 0.0, "upper": 1.0, "stoichiometry": {"a": 1.0}}})
    result = model.solve({}, set(), 5)
    assert type(result["numerical_pass"]) is bool
    assert json.loads(json.dumps(result, allow_nan=False))["numerical_pass"] is True


def test_cross_study_patient_is_fully_excluded():
    rows = [
        {"patient_id": "OCM74", "study_accession": "A"},
        {"patient_id": "OCM74", "study_accession": "B"},
        {"patient_id": "other", "study_accession": "B"},
        {"patient_id": "other", "study_accession": "B"},
    ]
    train, test = heldout_split(rows, "A", 3)
    assert test == [0] and len(train) == 1 and train[0] in (2, 3)
    assert one_per_patient(rows, range(4), 3) == one_per_patient(rows, range(4), 3)


def test_training_mad_does_not_use_heldout_extreme():
    train = np.array([[1, 2, 3], [1, 1, 1], [0, 1, 0]], dtype=float)
    assert select_mad(train, ["a", "b", "c"], 1).tolist() == [0]


def test_component_labels_are_scale_invariant():
    w = np.array([[1.0, 2.0], [3.0, 1.0]])
    h = np.array([[2.0, 1.0], [1.0, 2.0]])
    scale = np.array([100.0, 0.1])
    assert np.array_equal(normalized_labels(w, h), normalized_labels(w / scale, h * scale[:, None]))


def test_nested_gpr_and_missing_are_not_zero():
    rule = ("and", [("gene", "A"), ("or", [("gene", "B"), ("gene", "C")])])
    assert evaluate(rule, {"A": 3.0, "B": 1.0, "C": 1.0}) == 2.0
    assert evaluate(rule, {"A": 3.0, "B": 1.0}) is None
    assert evaluate(rule, {"A": 0.0, "B": 1.0, "C": 1.0}) == 0.0


def test_reverse_only_cap_and_no_rna_zero_knockout():
    reactions = {"r": {"lower": -10.0, "upper": 0.0}}
    caps, audit = cap_policy(reactions, {"r": ("gene", "A")}, {"A": 2.0}, 1.0)
    assert caps["r"] == (-2.0, 0.0)
    caps, audit = cap_policy(reactions, {"r": ("gene", "A")}, {"A": 0.0}, 1.0)
    assert not caps and audit["zero_rna_not_forced_to_knockout"] == 1


def test_uptake_direction_both_conventions():
    assert block_uptake({"stoichiometry": {"a": -1}, "lower": -10, "upper": 100}) == (0, 100)
    assert block_uptake({"stoichiometry": {"a": 1}, "lower": -100, "upper": 10}) == (-100, 0)
    with pytest.raises(ValueError):
        block_uptake({"stoichiometry": {"a": -1, "b": 1}, "lower": -10, "upper": 10})


def test_continuous_positive_control_and_boundary():
    reactions = {
        "uptake": {"lower": -10.0, "upper": 0.0, "stoichiometry": {"a": -1.0}},
        "enzyme": {"lower": 0.0, "upper": 100.0, "stoichiometry": {"a": -1.0, "b": 1.0}},
        "biomass_human": {"lower": 0.0, "upper": 100.0, "stoichiometry": {"b": -1.0}},
    }
    model = ContinuousModel(reactions)
    assert model.solve({}, set(), 5)["maximum_biomass"] == pytest.approx(10.0)
    assert model.solve({"enzyme": (0, 2)}, set(), 5)["maximum_biomass"] == pytest.approx(2.0)
    blocked = model.solve({}, {"uptake"}, 5)
    assert blocked["numerical_pass"] and blocked["maximum_biomass"] == pytest.approx(0.0)


def test_exact_signflip_and_bh():
    assert signflip([1, 1, 1])[0] == pytest.approx(0.25)
    assert bh([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.04, 0.04])


def test_nmf_training_convergence_on_small_positive_control():
    rng = np.random.default_rng(17)
    x = rng.uniform(0.1, 1, (12, 2)) @ rng.uniform(0.1, 1, (2, 20))
    model, w, diagnostics = fit_best(x, 2, 3, 17, 3000)
    assert any(d["converged"] for d in diagnostics)
    assert np.linalg.norm(x - w @ model.components_) / np.linalg.norm(x) < 0.02


def test_stroma_comparison_counts_patients_not_libraries():
    names = ["G" + str(i) for i in range(20)]
    x = np.tile(np.arange(20, dtype=float)[:, None], (1, 5))
    x[:10, [0, 1, 3]] += 30
    rows = [
        {
            "patient_id": patient,
            "study_accession": "S",
            "run_accession": str(i),
            "primary_cohort_eligible": "true" if kind == "tumour" else "false",
            "sample_class": kind,
        }
        for i, (patient, kind) in enumerate(
            [
                ("P1", "tumour"),
                ("P1", "tumour"),
                ("P1", "stroma"),
                ("P2", "tumour"),
                ("P2", "stroma"),
            ]
        )
    ]
    output = paired_stroma(x, rows, names, {"TEST": {"genes": set(names[:10])}})
    result = output["tests"][0]
    assert result["patient_count"] == 2
    assert result["matched_patient_study_count"] == 2
    assert result["mean_tumour_minus_stroma_percentile_score"] > 0
