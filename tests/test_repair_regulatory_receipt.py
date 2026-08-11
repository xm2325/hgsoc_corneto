from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "repair_regulatory_receipt.py"
SPEC = importlib.util.spec_from_file_location("repair_regulatory_receipt", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _root(samples):
    return {
        "schema_version": "corneto_regulatory_pilot.v1",
        "status": "partial",
        "study": "STUDY",
        "primary_only": True,
        "method": {"lambda_reg": 0.01, "response_blind": True},
        "source_sha256": {"expression": "abc"},
        "samples": samples,
    }


def _write(path: Path, value):
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_repair_replaces_only_errored_sample(tmp_path):
    original = _write(
        tmp_path / "original.json",
        _root(
            [
                {"run_accession": "A", "status": "optimal", "selected_edges": []},
                {"run_accession": "B", "status": "error", "selected_edges": []},
            ]
        ),
    )
    retry = _write(
        tmp_path / "retry.json",
        _root([{"run_accession": "B", "status": "optimal", "selected_edges": []}]),
    )
    result = MODULE.repair(original, [retry])
    assert result["status"] == "completed"
    assert result["schema_version"] == "corneto_regulatory_pilot_repaired.v1"
    assert [row["run_accession"] for row in result["samples"]] == ["A", "B"]
    assert result["samples"][1]["status"] == "optimal"
    assert result["repair"]["replacements"][0]["run_accession"] == "B"


def test_repair_rejects_non_error_replacement(tmp_path):
    original = _write(
        tmp_path / "original.json",
        _root([{"run_accession": "A", "status": "optimal", "selected_edges": []}]),
    )
    retry = _write(
        tmp_path / "retry.json",
        _root([{"run_accession": "A", "status": "optimal", "selected_edges": []}]),
    )
    with pytest.raises(MODULE.RepairError, match="is not errored"):
        MODULE.repair(original, [retry])


def test_repair_rejects_provenance_mismatch(tmp_path):
    original = _write(
        tmp_path / "original.json",
        _root([{"run_accession": "A", "status": "error", "selected_edges": []}]),
    )
    changed = _root([{"run_accession": "A", "status": "optimal", "selected_edges": []}])
    changed["source_sha256"] = {"expression": "changed"}
    retry = _write(tmp_path / "retry.json", changed)
    with pytest.raises(MODULE.RepairError, match="source_sha256 differs"):
        MODULE.repair(original, [retry])
