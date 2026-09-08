"""Strict I/O and statistical helpers for the post-audit diagnostic suite."""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

STUDIES = ("E-MTAB-7223", "E-MTAB-10801", "E-MTAB-11000", "E-MTAB-14568")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temp.replace(path)


def stamp():
    return datetime.now(UTC).isoformat()


def tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_data(root, rna_root):
    genes = names = None
    values, metadata, sources = [], [], []
    for study in STUDIES:
        matrix_path = Path(rna_root) / study / "gene_log1p_tpm.tsv.gz"
        qc_path = Path(root) / "data/processed/rna" / study / "aggregation/sample_qc.tsv"
        rows = tsv(qc_path)
        qc = {r["run_accession"]: r for r in rows}
        if len(qc) != len(rows):
            raise ValueError("duplicate QC run")
        with gzip.open(matrix_path, "rt") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            if header[:2] != ["gene_id", "gene_name"] or len(set(header[2:])) != len(header[2:]):
                raise ValueError("invalid expression header")
            records = list(reader)
        if any(len(row) != len(header) for row in records):
            raise ValueError("ragged expression matrix")
        ids = [r[0] for r in records]
        symbols = [r[1] for r in records]
        if len(set(ids)) != len(ids) or set(header[2:]) != set(qc):
            raise ValueError("gene/run identity mismatch")
        if genes is not None and (ids != genes or symbols != names):
            raise ValueError("cohort gene references differ")
        genes, names = ids, symbols
        x = np.asarray([r[2:] for r in records], dtype=float)
        if not np.isfinite(x).all() or (x < 0).any():
            raise ValueError("nonfinite/negative log1p TPM")
        for run in header[2:]:
            row = qc[run]
            if row["study_accession"] != study:
                raise ValueError("wrong study in QC")
            if row["primary_cohort_eligible"] == "true":
                if row["sample_class"] != "tumour" or row["histotype_group"] != "HGSOC":
                    raise ValueError("invalid primary membership")
                if row["patient_id"] in ("", "NA"):
                    raise ValueError("primary patient missing")
            metadata.append(row)
        values.append(x)
        sources.extend({"path": str(p), "sha256": sha(p)} for p in (matrix_path, qc_path))
    if len({m["run_accession"] for m in metadata}) != len(metadata):
        raise ValueError("duplicate cross-study run")
    primary = [m for m in metadata if m["primary_cohort_eligible"] == "true"]
    if len(primary) != 60 or len({m["patient_id"] for m in primary}) != 52:
        raise ValueError("frozen primary denominator must be 60 OCM / 52 patients")
    return genes, names, np.concatenate(values, axis=1), metadata, sources


def primary_data(x, metadata):
    keep = [i for i, m in enumerate(metadata) if m["primary_cohort_eligible"] == "true"]
    return x[:, keep], [metadata[i] for i in keep]


def select_mad(x, genes, count):
    mad = np.median(np.abs(x - np.median(x, axis=1)[:, None]), axis=1)
    eligible = [i for i in range(len(genes)) if mad[i] > 0]
    if len(eligible) < count:
        raise ValueError("insufficient positive-MAD training genes")
    return np.asarray(sorted(eligible, key=lambda i: (-mad[i], genes[i]))[:count])


def one_per_patient(metadata, indices, seed):
    rng = np.random.default_rng(seed)
    groups = {}
    for i in indices:
        groups.setdefault(metadata[i]["patient_id"], []).append(i)
    return sorted(int(rng.choice(sorted(groups[p]))) for p in sorted(groups))


def heldout_split(metadata, study, seed):
    test = [i for i, m in enumerate(metadata) if m["study_accession"] == study]
    patients = {metadata[i]["patient_id"] for i in test}
    candidates = [i for i, m in enumerate(metadata) if m["patient_id"] not in patients]
    train = one_per_patient(metadata, candidates, seed)
    if not train or not test:
        raise ValueError("empty grouped split")
    assert not ({metadata[i]["patient_id"] for i in train} & patients)
    return train, test


def bh(pvalues):
    values = np.asarray(pvalues, dtype=float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("invalid p-values")
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    result = np.empty(len(values))
    result[order] = adjusted
    return result.tolist()


def signflip(values, permutations=9999, seed=20260908):
    """Two-sided family-level test; exact for <=16 independent families."""
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("invalid sign-flip values")
    observed = abs(float(values.mean()))
    if len(values) <= 16:
        signs = np.asarray(list(itertools.product((-1, 1), repeat=len(values))))
        exceed = np.sum(np.abs(signs @ values / len(values)) >= observed - 1e-12)
        return float(exceed / len(signs)), "exact_family_signflip"
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1, 1], size=(permutations, len(values)))
    exceed = np.sum(np.abs(signs @ values / len(values)) >= observed - 1e-12)
    return float((exceed + 1) / (permutations + 1)), "monte_carlo_family_signflip"
