#!/usr/bin/env python3
"""Exploratory biological annotation with explicit denominators and controls.

Hallmark over-representation annotates NMF loadings, not independent DE.
Rank-based scores are not GSVA/GSEA. Patient/stroma comparisons and within-
patient distances are exploratory; no clinical response or causality is used.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from research_validation_common import (
    bh,
    load_data,
    primary_data,
    select_mad,
    sha,
    signflip,
    stamp,
    write_json,
)
from scipy.spatial.distance import pdist, squareform
from scipy.stats import hypergeom, rankdata


def read_gmt(path):
    sets = {}
    with Path(path).open() as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < 3 or not row[0].startswith("HALLMARK_") or row[0] in sets:
                raise ValueError("invalid Hallmark GMT")
            sets[row[0]] = {"url": row[1], "genes": set(row[2:])}
    if len(sets) != 50:
        raise ValueError("expected the 50-set Hallmark collection")
    return sets


def loading_annotation(root, symbol_map, sets, top_count=200):
    tests, components, sources = [], [], []
    for rank in (2, 3):
        path = (
            root
            / "data/processed/rna/pooled_primary60/nmf"
            / f"rank_{rank}/best_w_gene_loadings.tsv.gz"
        )
        sources.append({"path": str(path), "sha256": sha(path)})
        with gzip.open(path, "rt") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            rows = list(reader)
        if len(header) != rank + 1:
            raise ValueError("loading dimension mismatch")
        if any(r[0] not in symbol_map for r in rows):
            raise ValueError("loading gene absent from expression reference")
        weights = np.asarray([r[1:] for r in rows], dtype=float)
        if (
            not np.isfinite(weights).all()
            or (weights < 0).any()
            or (weights.sum(axis=0) <= 0).any()
        ):
            raise ValueError("invalid NMF weights")
        weights = weights / weights.sum(axis=0)[None, :]
        background = {symbol_map[r[0]] for r in rows if symbol_map[r[0]] not in ("", "NA")}
        for component in range(rank):
            specificity = weights[:, component] - np.mean(
                np.delete(weights, component, axis=1), axis=1
            )
            order = sorted(range(len(rows)), key=lambda i: (-specificity[i], rows[i][0]))
            selected, top_genes = set(), []
            for i in order:
                symbol = symbol_map[rows[i][0]]
                if specificity[i] <= 0:
                    break
                if symbol not in background or symbol in selected:
                    continue
                selected.add(symbol)
                top_genes.append(
                    {
                        "gene_id": rows[i][0],
                        "symbol": symbol,
                        "loading_specificity": float(specificity[i]),
                    }
                )
                if len(selected) == top_count:
                    break
            if len(selected) < 20:
                raise ValueError("insufficient component-specific genes")
            components.append(
                {
                    "rank": rank,
                    "state": header[component + 1],
                    "top_genes": top_genes,
                    "background_symbol_count": len(background),
                    "selected_symbol_count": len(selected),
                }
            )
            for name, gene_set in sets.items():
                background_set = gene_set["genes"] & background
                hits = sorted(selected & background_set)
                p = float(
                    hypergeom.sf(len(hits) - 1, len(background), len(background_set), len(selected))
                )
                tests.append(
                    {
                        "rank": rank,
                        "state": header[component + 1],
                        "hallmark": name,
                        "background_n": len(background),
                        "set_in_background_n": len(background_set),
                        "selected_n": len(selected),
                        "overlap_n": len(hits),
                        "overlap_genes": hits,
                        "p_value": p,
                        "gene_set_url": gene_set["url"],
                    }
                )
    for row, q in zip(tests, bh([r["p_value"] for r in tests]), strict=True):
        row["q_value_bh_all_rank_state_sets"] = q
    return components, tests, sources


def family_distance(x, metadata, genes):
    x, meta = primary_data(x, metadata)
    chosen = select_mad(x, genes, min(6000, len(genes)))
    distances = squareform(pdist(x[chosen].T, metric="correlation"))
    if not np.isfinite(distances).all():
        raise ValueError("undefined expression correlation distance")
    families = defaultdict(list)
    for i, row in enumerate(meta):
        families[row["patient_id"]].append(i)
    all_pairs = list(itertools.combinations(range(len(meta)), 2))
    records = []
    for patient, indices in sorted(families.items()):
        if len(indices) < 2:
            continue
        within = list(itertools.combinations(indices, 2))
        matched_means, comparisons = [], []
        for a, b in within:
            study_pair = sorted([meta[a]["study_accession"], meta[b]["study_accession"]])
            eligible = [
                (i, j)
                for i, j in all_pairs
                if meta[i]["patient_id"] != meta[j]["patient_id"]
                and patient not in (meta[i]["patient_id"], meta[j]["patient_id"])
                and sorted([meta[i]["study_accession"], meta[j]["study_accession"]]) == study_pair
            ]
            if not eligible:
                raise ValueError("no matched between-patient reference pairs")
            # Average first within each reference patient pair, preventing
            # patients with multiple OCMs from dominating this baseline.
            by_pair = defaultdict(list)
            for i, j in eligible:
                by_pair[tuple(sorted([meta[i]["patient_id"], meta[j]["patient_id"]]))].append(
                    distances[i, j]
                )
            mean = float(np.mean([np.mean(v) for v in by_pair.values()]))
            matched_means.append(mean)
            comparisons.append(
                {
                    "runs": [meta[a]["run_accession"], meta[b]["run_accession"]],
                    "within_distance": float(distances[a, b]),
                    "matched_between_distance": mean,
                    "matched_patient_pair_count": len(by_pair),
                }
            )
        observed = float(np.mean([distances[a, b] for a, b in within]))
        reference = float(np.mean(matched_means))
        records.append(
            {
                "patient_id": patient,
                "ocms": [meta[i]["canonical_ocm_id"] for i in indices],
                "comparisons": comparisons,
                "mean_within_distance": observed,
                "mean_matched_between_distance": reference,
                "within_minus_between": observed - reference,
            }
        )
    # Shared between-patient reference sets induce dependence: report no
    # naive pair-level p-value; this is descriptive family-level evidence.
    return {
        "families": records,
        "repeated_patient_count": len(records),
        "claim_limit": (
            "Pooled-MAD descriptive distances; shared reference pairs and only seven families "
            "preclude treating OCM pairs as independent replicates. "
            "No temporal or acquired-resistance claim."
        ),
    }


def paired_stroma(x, metadata, names, sets):
    # Duplicate symbols are averaged before within-sample percentile ranking.
    groups = defaultdict(list)
    for i, name in enumerate(names):
        if name not in ("", "NA"):
            groups[name].append(i)
    symbols = sorted(groups)
    matrix = np.asarray([x[groups[s]].mean(axis=0) for s in symbols])
    ranks = rankdata(matrix, axis=0, method="average") / len(symbols)
    symbol_index = {s: i for i, s in enumerate(symbols)}
    pairs = defaultdict(lambda: {"tumour": [], "stroma": []})
    for i, row in enumerate(metadata):
        if row.get("patient_id") in (None, "", "NA"):
            continue
        if row["primary_cohort_eligible"] == "true":
            pairs[(row["patient_id"], row["study_accession"])]["tumour"].append(i)
        elif row["sample_class"] == "stroma":
            pairs[(row["patient_id"], row["study_accession"])]["stroma"].append(i)
    complete = {k: v for k, v in pairs.items() if v["tumour"] and v["stroma"]}
    tests = []
    for name, gene_set in sets.items():
        keep = [symbol_index[g] for g in sorted(gene_set["genes"] & set(symbol_index))]
        if len(keep) < 10:
            continue
        score = ranks[keep].mean(axis=0)
        per_patient = defaultdict(list)
        for (patient, _study), indices in complete.items():
            per_patient[patient].append(
                float(score[indices["tumour"]].mean() - score[indices["stroma"]].mean())
            )
        values = [float(np.mean(v)) for patient, v in sorted(per_patient.items())]
        if not values:
            continue
        p, method = signflip(values)
        rng = np.random.default_rng(20260908)
        boot = np.asarray(values)[rng.integers(0, len(values), size=(2000, len(values)))].mean(
            axis=1
        )
        tests.append(
            {
                "hallmark": name,
                "gene_coverage": len(keep),
                "set_size": len(gene_set["genes"]),
                "patient_count": len(values),
                "matched_patient_study_count": len(complete),
                "mean_tumour_minus_stroma_percentile_score": float(np.mean(values)),
                "patient_bootstrap_95_interval": np.quantile(boot, [0.025, 0.975]).tolist(),
                "patient_differences": [
                    {"patient_id": p, "difference": float(np.mean(v))}
                    for p, v in sorted(per_patient.items())
                ],
                "p_value": p,
                "test": method,
            }
        )
    for row, q in zip(tests, bh([r["p_value"] for r in tests]), strict=True):
        row["q_value_bh"] = q
    return {
        "tests": tests,
        "matched_pairs": [
            {
                "patient_id": k[0],
                "study": k[1],
                "tumour_runs": [metadata[i]["run_accession"] for i in v["tumour"]],
                "stroma_runs": [metadata[i]["run_accession"] for i in v["stroma"]],
            }
            for k, v in sorted(complete.items())
        ],
        "claim_limit": (
            "Exploratory same-study patient-matched RNA programme difference; mean percentile-rank "
            "score is not pathway flux, GSVA or GSEA. "
            "Lineage/proliferation/culture effects remain plausible explanations."
        ),
    }


def regulatory_description(root):
    base = root / "data/processed/corneto/regulatory_multisample/v1"
    receipt_path, bundle_path = base / "grid/pooled/l0p001.json", base / "input_bundle.json"
    receipt, bundle = [json.loads(p.read_text()) for p in (receipt_path, bundle_path)]
    if receipt.get("status") != "completed" or len(receipt["conditions"]) != 60:
        raise ValueError("regulatory reference does not pass expected contract")
    counts, patients = Counter(), defaultdict(set)
    for condition in receipt["conditions"]:
        if condition["status"] != "optimal":
            raise ValueError("nonoptimal regulatory reference condition")
        for edge in condition["selected_edges"]:
            key = edge["source"], edge["target"], edge["sign"]
            counts[key] += 1
            patients[key].add(condition["patient_id"])
    targets = Counter(tf for c in bundle["conditions"] for tf in c["outputs"])
    inputs = Counter(tf for c in bundle["conditions"] for tf in c["inputs"])
    return {
        "selected_edges": [
            {
                "source": k[0],
                "target": k[1],
                "sign": k[2],
                "ocm_count": n,
                "distinct_patient_count": len(patients[k]),
            }
            for k, n in counts.most_common()
        ],
        "prespecified_output_frequency": dict(targets.most_common()),
        "input_frequency": dict(inputs.most_common()),
        "sources": [{"path": str(p), "sha256": sha(p)} for p in (receipt_path, bundle_path)],
        "claim_limit": (
            "Selected signed edges under the existing pooled prior/constraints; output-TF "
            "prominence is partly by design. Not binding, causality, phosphorylation "
            "or external replication."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--rna-root", type=Path, required=True)
    parser.add_argument("--hallmark-gmt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "existing_biology_interpretation.v1",
        "status": "running",
        "started_at_utc": stamp(),
        "script_sha256": sha(__file__),
        "common_code_sha256": sha(Path(__file__).with_name("research_validation_common.py")),
        "gurobi_sessions_requested": 0,
        "hallmark_source": {
            "url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt",
            "version": "2025.1.Hs",
            "sha256": sha(args.hallmark_gmt),
        },
        "claim_limit": (
            "Exploratory biological annotation of existing expression/regulatory results. "
            "No metabolic partial incumbent is interpreted biologically."
        ),
    }
    output = args.output_dir / "receipt.json"
    write_json(output, report)
    try:
        genes, names, x, metadata, sources = load_data(args.root, args.rna_root)
        report["sources"] = sources
        sets = read_gmt(args.hallmark_gmt)
        components, enrichment, sources = loading_annotation(
            args.root, dict(zip(genes, names, strict=True)), sets
        )
        report.update(nmf_components=components, nmf_hallmark_annotation=enrichment)
        report["sources"].extend(sources)
        report["within_patient_expression"] = family_distance(x, metadata, genes)
        report["matched_tumour_stroma"] = paired_stroma(x, metadata, names, sets)
        report["regulatory"] = regulatory_description(args.root)
        report["status"] = "completed"
    except Exception as exc:
        report["status"], report["error"] = "incomplete", repr(exc)
        raise
    finally:
        report["finished_at_utc"] = stamp()
        write_json(output, report)
    print("completed", output, flush=True)


if __name__ == "__main__":
    main()
