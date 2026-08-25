from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", dtype=str, comment=None)
    df.columns = [c.lstrip("#") for c in df.columns]
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
