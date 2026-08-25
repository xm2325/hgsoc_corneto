from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    """Read a PLINK text table while preserving its real header.

    PLINK .pvar files may contain VCF-style ``##`` metadata records before the
    tabular ``#CHROM`` header. Other PLINK reports normally begin directly
    with a single-hash header such as ``#FID``. We therefore locate the first
    single-hash header explicitly instead of letting pandas infer a header from
    the first physical line.
    """
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    header_idx = next(
        (idx for idx, line in enumerate(lines) if line.startswith("#") and not line.startswith("##")),
        None,
    )
    if header_idx is None:
        raise ValueError(f"no PLINK table header found in {path}")

    df = pd.read_csv(io.StringIO("".join(lines[header_idx:])), sep=r"\s+", dtype=str)
    df.columns = [column.lstrip("#") for column in df.columns]
    return df


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        raise ValueError(f"refusing to write empty table: {path.name}")
    df.to_parquet(path, index=False)


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
    for name, source in sources.items():
        frame = read_table(source)
        target = outdir / f"{name}.parquet"
        write_parquet(frame, target)
        counts[name] = int(len(frame))
        files[name] = target.name

    if counts["variants"] <= 0 or counts["samples"] <= 0:
        raise ValueError("QC output must contain at least one variant and one sample")

    summary = {
        "sample_count": counts["samples"],
        "variant_count": counts["variants"],
        "row_counts": counts,
        "parquet_files": files,
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    for name in ("pvar", "psam", "afreq", "vmiss", "smiss", "hardy", "eigenvec"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--summary", required=True)
    return p


if __name__ == "__main__":
    build_outputs(parser().parse_args())
