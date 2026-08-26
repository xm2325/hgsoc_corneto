from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import numbers
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from pandas.api.types import is_bool_dtype, is_float_dtype, is_integer_dtype


INTEGER_COLUMNS = {
    "POS",
    "OBS_CT",
    "MISSING_CT",
    "HOM_A1_CT",
    "HET_A1_CT",
    "TWO_AX_CT",
}
FLOAT_COLUMNS = {"ALT_FREQS", "F_MISS", "O_HET", "E_HET", "P"}


def read_table(path: Path) -> pd.DataFrame:
    """Read a PLINK text table while preserving its real header."""
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    header_idx = next(
        (idx for idx, line in enumerate(lines) if line.startswith("#") and not line.startswith("##")),
        None,
    )
    if header_idx is None:
        raise ValueError(f"no PLINK table header found in {path}")

    df = pd.read_csv(
        io.StringIO("".join(lines[header_idx:])),
        sep=r"\s+",
        dtype=str,
        keep_default_na=False,
    )
    df.columns = [column.lstrip("#") for column in df.columns]
    return df


def coerce_types(table_name: str, df: pd.DataFrame) -> pd.DataFrame:
    """Convert PLINK text columns into stable analysis-layer data types."""
    typed = df.copy()
    for column in typed.columns:
        if column in INTEGER_COLUMNS:
            typed[column] = pd.to_numeric(typed[column], errors="raise").astype("int64")
        elif column in FLOAT_COLUMNS or column.startswith("PC"):
            typed[column] = pd.to_numeric(typed[column], errors="raise").astype("float64")
        else:
            typed[column] = typed[column].astype("string")

    if table_name == "variants" and {"CHROM", "POS", "REF", "ALT"}.issubset(typed.columns):
        typed = typed.sort_values(["CHROM", "POS", "REF", "ALT"], kind="stable").reset_index(drop=True)
    return typed


def _semantic_type(series: pd.Series) -> str:
    if is_bool_dtype(series.dtype):
        return "boolean"
    if is_integer_dtype(series.dtype):
        return "integer"
    if is_float_dtype(series.dtype):
        return "float"
    return "string"


def _canonical_scalar(value: object) -> list[object]:
    if value is None or pd.isna(value):
        return ["null", None]
    if isinstance(value, bool):
        return ["boolean", bool(value)]
    if isinstance(value, numbers.Integral):
        return ["integer", str(int(value))]
    if isinstance(value, numbers.Real):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError("semantic hashing rejects non-finite floating-point values")
        if number == 0.0:
            number = 0.0
        return ["float", number.hex()]
    return ["string", str(value)]


def semantic_table_sha256(df: pd.DataFrame) -> str:
    """Hash logical typed table content independently of Parquet encoding bytes.

    Row and column order are part of the contract. This is intentional: sample order
    and the pipeline's deterministic table ordering are release semantics, while
    compression, row groups and file metadata are not.
    """
    columns = [
        {"name": str(column), "semantic_type": _semantic_type(df[column])}
        for column in df.columns
    ]
    digest = hashlib.sha256()
    digest.update(
        json.dumps({"columns": columns}, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )
    for row in df.itertuples(index=False, name=None):
        canonical_row = [_canonical_scalar(value) for value in row]
        digest.update(b"\n")
        digest.update(
            json.dumps(canonical_row, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        )
    return digest.hexdigest()


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        raise ValueError(f"refusing to write empty table: {path.name}")
    df.to_parquet(path, index=False, compression="zstd")


def parquet_schema(path: Path) -> dict[str, object]:
    parquet_file = pq.ParquetFile(path)
    schema = parquet_file.schema_arrow
    return {
        "row_count": int(parquet_file.metadata.num_rows),
        "columns": [
            {"name": field.name, "type": str(field.type), "nullable": field.nullable}
            for field in schema
        ],
    }


def build_outputs(args: argparse.Namespace) -> dict[str, object]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sources = {
        "variants": Path(args.pvar),
        "samples": Path(args.psam),
        "allele_frequencies": Path(args.afreq),
        "variant_missingness": Path(args.vmiss),
        "sample_missingness": Path(args.smiss),
        "hardy_weinberg": Path(args.hardy),
        "pca_scores": Path(args.eigenvec),
    }

    counts: dict[str, int] = {}
    files: dict[str, str] = {}
    schemas: dict[str, object] = {}
    semantic_hashes: dict[str, str] = {}
    for name, source in sources.items():
        frame = coerce_types(name, read_table(source))
        target = outdir / f"{name}.parquet"
        write_parquet(frame, target)
        counts[name] = int(len(frame))
        files[name] = target.name
        schemas[name] = parquet_schema(target)
        semantic_hashes[name] = semantic_table_sha256(frame)

    if counts["variants"] <= 0 or counts["samples"] <= 0:
        raise ValueError("QC output must contain at least one variant and one sample")

    summary = {
        "sample_count": counts["samples"],
        "variant_count": counts["variants"],
        "row_counts": counts,
        "parquet_files": files,
        "semantic_hashes": semantic_hashes,
        "storage": {"format": "parquet", "compression": "zstd"},
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    Path(args.schema_manifest).write_text(
        json.dumps({"format": "parquet", "tables": schemas}, indent=2, sort_keys=True) + "\n"
    )
    return summary


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    for name in ("pvar", "psam", "afreq", "vmiss", "smiss", "hardy", "eigenvec"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--schema-manifest", required=True)
    return p


if __name__ == "__main__":
    build_outputs(parser().parse_args())
