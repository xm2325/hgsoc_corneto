import json
import sys

import numpy as np
import pytest
import validate_binary_support as runner


@pytest.mark.parametrize("native", [False, True])
def test_binary_runner_and_independent_support(tmp_path, monkeypatch, native):
    reactions = {
        "in": {"lower": 0, "upper": 10, "stoichiometry": {"a": 1}},
        "biomass_human": {"lower": 0, "upper": 10, "stoichiometry": {"a": -1}},
    }
    model = tmp_path / "model"
    model.write_text("test")
    source = tmp_path / "data/processed/revised_gpr_pfba_20260910/receipt.json"
    source.parent.mkdir(parents=True)
    source.write_text(
        json.dumps(
            {
                "status": "completed",
                "model_sha256": runner.sha(model),
                "cases": [
                    {
                        "fraction": 0.9,
                        "condition": "test",
                        "selected_reactions": list(reactions),
                        "growth_maximum": {"maximum_biomass": 10},
                        "flux": {"in": 9.0, "biomass_human": 9.0},
                    }
                ],
            }
        )
    )
    monkeypatch.setattr(runner, "read_sbml", lambda p: reactions)
    monkeypatch.setattr(runner, "read_gprs", lambda p: ({}, []))
    monkeypatch.setattr(
        runner,
        "load_data",
        lambda *a: (["ENSG1.1"], [], np.array([[1.0]]), [{"run_accession": "test"}], []),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "test",
            "--root",
            str(tmp_path),
            "--rna",
            str(tmp_path),
            "--model",
            str(model),
            "--output",
            str(tmp_path / "out"),
        ],
    )
    if native:
        pytest.importorskip("gurobipy")
        sys.argv.append("--native-indicators")
    runner.main()
    d = json.loads((tmp_path / "out/receipt.json").read_text())
    assert d["status"] == "validated_restricted_optimum"
    assert d["unselected_flux"] == 0
    assert d["fixed_support_check"]["maximum_biomass"] == 10
