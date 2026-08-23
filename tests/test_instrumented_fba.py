from types import SimpleNamespace

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
