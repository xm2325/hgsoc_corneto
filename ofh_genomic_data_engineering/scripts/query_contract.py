from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd


def _column_types(schema_manifest: dict[str, object], table: str) -> dict[str, str]:
    table_schema = schema_manifest["tables"][table]
    return {column["name"]: column["type"] for column in table_schema["columns"]}


def _write_payload(
    *,
    checks: list[dict[str, object]],
    output_path: Path,
    total_variants: int | None,
    region_query: dict[str, int] | None,
) -> dict[str, object]:
    passed = all(check["status"] == "PASS" for check in checks)
    payload = {
        "status": "PASS" if passed else "FAIL",
        "engine": {"name": "duckdb", "version": duckdb.__version__},
        "total_variants": total_variants,
        "region_query": region_query,
        "checks": checks,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def validate_query_layer(
    *, variants_path: Path, schema_manifest_path: Path, output_path: Path
) -> dict[str, object]:
    schema = json.loads(schema_manifest_path.read_text())
    checks: list[dict[str, object]] = []

    def record(name: str, passed: bool, detail: object) -> None:
        checks.append({"name": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    types = _column_types(schema, "variants")
    pos_ok = types.get("POS") == "int64"
    chrom_ok = types.get("CHROM") == "string"
    record("variants.POS_int64", pos_ok, types.get("POS"))
    record("variants.CHROM_string", chrom_ok, types.get("CHROM"))

    # Schema validation is a hard precondition for the SQL contract. Returning
    # a structured FAIL here avoids replacing a data-contract result with a
    # lower-level DuckDB binder exception when POS is incorrectly stored as text.
    if not (pos_ok and chrom_ok):
        return _write_payload(
            checks=checks,
            output_path=output_path,
            total_variants=None,
            region_query=None,
        )

    con = duckdb.connect(database=":memory:")
    stats = con.execute(
        "SELECT COUNT(*), MIN(POS), MAX(POS) FROM read_parquet(?)", [str(variants_path)]
    ).fetchone()
    total, min_pos, max_pos = int(stats[0]), int(stats[1]), int(stats[2])
    record("query.total_positive", total > 0, total)
    record("query.position_order", min_pos <= max_pos, {"min": min_pos, "max": max_pos})

    region_start = min_pos
    region_end = min(max_pos, min_pos + 1_000_000)
    duck_count = int(
        con.execute(
            "SELECT COUNT(*) FROM read_parquet(?) WHERE POS BETWEEN ? AND ?",
            [str(variants_path), region_start, region_end],
        ).fetchone()[0]
    )
    frame = pd.read_parquet(variants_path, columns=["POS"])
    pandas_count = int(frame["POS"].between(region_start, region_end, inclusive="both").sum())
    record(
        "query.region_count_matches",
        duck_count == pandas_count and duck_count > 0,
        {
            "start": region_start,
            "end": region_end,
            "duckdb": duck_count,
            "pandas": pandas_count,
        },
    )

    positions = [
        int(row[0])
        for row in con.execute(
            "SELECT POS FROM read_parquet(?) WHERE POS BETWEEN ? AND ? ORDER BY POS LIMIT 100",
            [str(variants_path), region_start, region_end],
        ).fetchall()
    ]
    record("query.ordered_positions", positions == sorted(positions), {"rows_checked": len(positions)})

    return _write_payload(
        checks=checks,
        output_path=output_path,
        total_variants=total,
        region_query={
            "start": region_start,
            "end": region_end,
            "matched_variants": duck_count,
        },
    )


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--variants", required=True)
    p.add_argument("--schema-manifest", required=True)
    p.add_argument("--output", required=True)
    return p


def main() -> int:
    args = parser().parse_args()
    payload = validate_query_layer(
        variants_path=Path(args.variants),
        schema_manifest_path=Path(args.schema_manifest),
        output_path=Path(args.output),
    )
    if payload["status"] != "PASS":
        failed = [check["name"] for check in payload["checks"] if check["status"] == "FAIL"]
        print("query contract failed: " + ", ".join(failed))
        return 2
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
