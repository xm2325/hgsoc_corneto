import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.query_contract import validate_query_layer


def _write_fixture(tmp_path: Path, pos_type: str = "int64") -> tuple[Path, Path]:
    pytest.importorskip("pyarrow")
    variants = pd.DataFrame(
        {
            "CHROM": pd.Series(["22", "22", "22"], dtype="string"),
            "POS": pd.Series([100, 200, 2_000_000], dtype=pos_type),
            "ID": pd.Series(["rs1", "rs2", "rs3"], dtype="string"),
            "REF": pd.Series(["A", "C", "G"], dtype="string"),
            "ALT": pd.Series(["G", "T", "A"], dtype="string"),
        }
    )
    variants_path = tmp_path / "variants.parquet"
    variants.to_parquet(variants_path, index=False, compression="zstd")
    schema_manifest = {
        "format": "parquet",
        "tables": {
            "variants": {
                "row_count": 3,
                "columns": [
                    {"name": "CHROM", "type": "string", "nullable": True},
                    {"name": "POS", "type": "int64" if pos_type == "int64" else "string", "nullable": True},
                    {"name": "ID", "type": "string", "nullable": True},
                    {"name": "REF", "type": "string", "nullable": True},
                    {"name": "ALT", "type": "string", "nullable": True},
                ],
            }
        },
    }
    schema_path = tmp_path / "schema_manifest.json"
    schema_path.write_text(json.dumps(schema_manifest))
    return variants_path, schema_path


def test_query_contract_passes_for_typed_parquet(tmp_path: Path) -> None:
    variants, schema = _write_fixture(tmp_path)
    payload = validate_query_layer(
        variants_path=variants,
        schema_manifest_path=schema,
        output_path=tmp_path / "query_validation.json",
    )
    assert payload["status"] == "PASS"
    assert payload["total_variants"] == 3
    assert payload["region_query"]["matched_variants"] == 2


def test_query_contract_rejects_string_position_schema(tmp_path: Path) -> None:
    variants, schema = _write_fixture(tmp_path, pos_type="string")
    payload = validate_query_layer(
        variants_path=variants,
        schema_manifest_path=schema,
        output_path=tmp_path / "query_validation.json",
    )
    assert payload["status"] == "FAIL"
    failed = {check["name"] for check in payload["checks"] if check["status"] == "FAIL"}
    assert "variants.POS_int64" in failed
