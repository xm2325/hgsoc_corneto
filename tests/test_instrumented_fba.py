from types import SimpleNamespace

import pytest

from hgsoc_corneto.metabolic.instrumented_fba import _status_name, _telemetry


class FakeModel:
    Status = 9
    SolCount = 1
    ObjVal = 12.0
    ObjBound = 10.0
    MIPGap = 1.0 / 6.0
    NodeCount = 42.0
    Runtime = 60.0
    Work = 4.0
    IterCount = 5.0
    BarIterCount = 6
    NumVars = 100
    NumBinVars = 20
    NumIntVars = 20
    NumConstrs = 80
    NumNZs = 400


def test_extracts_time_limited_incumbent_telemetry():
    solution = SimpleNamespace(
        status="user_limit",
        solver_stats=SimpleNamespace(extra_stats=FakeModel()),
    )
    telemetry, model = _telemetry(solution, max_seconds=60, mip_gap=1e-4, threads=8)
    assert model is not None
    assert telemetry.gurobi_status_name == "TIME_LIMIT"
    assert telemetry.has_incumbent is True
    assert telemetry.objective_value == 12.0
    assert telemetry.best_bound == 10.0
    assert telemetry.absolute_gap == 2.0
    assert telemetry.relative_gap == 1.0 / 6.0
    assert telemetry.variable_count == 100


def test_maps_known_and_unknown_status_codes():
    assert _status_name(9) == "TIME_LIMIT"
    assert _status_name(999) == "STATUS_999"


@pytest.mark.parametrize("indicators,expected", [((1.0, 1.0), True), ((9e-6, 9e-6), False)])
def test_finish_never_promotes_numerically_invalid_selection(tmp_path, indicators, expected):
    from test_metabolic_validation import primal_model

    from hgsoc_corneto.metabolic.instrumented_fba import _finish

    class CompleteModel(FakeModel):
        Status = 2
        ObjVal = -10.0
        ObjBound = -10.0
        MIPGap = 0.0

        def write(self, filename):
            from pathlib import Path

            Path(filename).write_text("mock incumbent artifact\n")

    model, problem = primal_model((0.009, 0.009), indicators)
    solution = SimpleNamespace(
        status="optimal",
        value=-10.0,
        solver_stats=SimpleNamespace(extra_stats=CompleteModel()),
    )
    log = tmp_path / "solver.log"
    log.write_text("mock optimal solver log\n")
    result = _finish(
        problem=problem,
        solution=solution,
        reaction_ids=["in", "out"],
        conditions=["A"],
        max_seconds=60,
        mip_gap=1e-4,
        threads=8,
        active_tolerance=1e-7,
        artifact_prefix=tmp_path / "solution",
        log_file=log,
        cobra_model=model,
        reaction_bounds={"A": {}},
    )
    assert result.scientific_success is expected
    assert result.status == ("completed" if expected else "partial_incumbent")
    assert len(result.artifact_sha256) == 3
