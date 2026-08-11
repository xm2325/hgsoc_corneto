# HGSOC × CORNETO

Reproducible analysis of patient-derived ovarian cancer models (OCMs), with two
separate scientific tracks:

1. label-blind regulatory-network inference followed by paclitaxel-response
   association testing; and
2. a methodological comparison of sequential, global, and joint metabolic
   constraint-based modelling.

The analysis contract is frozen in [`config/analysis_contract.yaml`](config/analysis_contract.yaml).
The primary cohort is **HGSOC tumour-only**. Paclitaxel AUC is the primary
outcome and log10(GI50) is secondary. Cross-sectional, intrinsic, and acquired
resistance are separate estimands, and every predictive evaluation is split by
patient.

## Current reproducible assets

- BioStudies IDF/SDRF snapshots for `E-MTAB-7223`, `E-MTAB-10801`,
  `E-MTAB-11000`, and `E-MTAB-14568`.
- ENA run-level reports for 117 paired-end RNA-seq runs (234 FASTQ files;
  444.951 GiB compressed).
- A parser for Tighe et al. Table S1 (83 OCMs from 68 patients), including
  histotype, biopsy type, longitudinal family, and chemotherapy-naive status.
- A parser for the public ABCB1 Table S2 workbook.
- A normalized run → library → OCM → patient manifest with explicit provenance
  and non-guessing alias rules.

The exact OCM-level paclitaxel AUC/GI50 and cumulative exposure values are not
contained in the public supplementary spreadsheets. They remain intentionally
missing rather than being digitized from figures.

## Quick start

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev,metadata]'
.venv/bin/python scripts/build_manifest.py \
  --tighe-pdf tmp/pdfs/mmc1.pdf
.venv/bin/python -m pytest
```

The RNA workflow is designed for Slurm/CSF3 because the public FASTQs do not fit
on a typical laptop. See [`docs/rna_reproduction.md`](docs/rna_reproduction.md)
before starting any transfer.

## Scientific boundaries

- Expression is observed; TF activity is inferred.
- A CORNETO network is a model-supported hypothesis, not a Bayesian posterior.
- Drug response is never supplied to network inference.
- Treatment-associated rewiring is not called causal without perturbation.
- Single edges are not treated as stable discoveries; modules and solution
  ensembles are the reporting unit.

## Sources

The source registry and immutable upstream revisions are recorded in
[`config/sources.yaml`](config/sources.yaml). Derived tables carry source fields,
and raw metadata snapshots are checksummed.
