# Taylor/HGSOC three-track analysis preparation

This repository freezes the analysis design for the Taylor-group HGSOC RNA
cohort without inventing or imputing paclitaxel phenotypes.  The preparation
script is metadata-only: it validates the cohort manifest, identifies the
longitudinal families that can support the acquired track, and assigns
deterministic patient-grouped cross-validation folds.

Run from the repository root:

```bash
python scripts/prepare_taylor_analysis.py \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --families data/processed/metadata/longitudinal_families.tsv \
  --schema config/taylor_analysis.json \
  --output-dir data/processed/taylor
```

Use `--overwrite` only when intentionally replacing a previous preparation
directory.  The command reads the two metadata TSVs and the JSON schema; it
does not require the phenotype tables and does not populate their columns.

## Frozen cohort and outcomes

The manifest is expected to contain 117 RNA-run rows.  The primary cohort is
the 60-row tumour-only HGSOC set with one representative RNA library per OCM,
covering 60 OCMs and 52 patients.  The current study distribution in that
primary set is E-MTAB-7223 (9), E-MTAB-10801 (13), E-MTAB-11000 (11), and
E-MTAB-14568 (27).  These are validation targets, not a replacement for the
manifest checks.

Paclitaxel AUC is the primary outcome and exact GI50 in nM is the secondary
outcome, analysed as `log10(GI50)`.  An outcome is usable only when its exact
availability flag is `true` and the corresponding numeric value is present and
valid.  A populated value with a `false` flag is rejected.  Consequently, the
preparation outputs currently carry availability flags only; they do not claim
that a Tighe phenotype value is available.

## Three tracks

### Cross-sectional

The unit is the representative tumour OCM RNA library in the 60-row primary
cohort.  It estimates baseline feature associations with contemporaneous
paclitaxel response.  Patient IDs, rather than OCM IDs, define the five CV
groups, so repeated OCMs from one patient cannot cross a train/test boundary.

### Intrinsic

This track is restricted to primary-cohort rows with an explicitly measured
cumulative paclitaxel exposure equal to zero.  `chemo_naive_at_biopsy` is not
used as a proxy for exact zero exposure.  Until the exact exposure and response
tables are available, the track remains blocked and the script emits no
intrinsic fold rows.

### Acquired

This track uses only `relationship_type == longitudinal` families in
`longitudinal_families.tsv`.  The family design output distinguishes:

- `pair_ready`: at least two family OCMs are in the primary cohort;
- `partial_primary_family`: one family OCM is in the primary cohort;
- `not_in_primary_cohort`: no family OCM is in the primary cohort.

The current metadata contain seven pair-ready primary longitudinal families,
two partial primary anchors (OCM231 and OCM341), and two longitudinal families
absent from the primary cohort (OCM118 and OCM124).  Mixed/spatial relationship
rows are not silently treated as longitudinal.  The acquired estimand is
within-patient change in features versus change in paclitaxel response.  Any
interpretation is limited to treatment-associated rewiring; the design does
not establish that treatment caused the change.

## Outputs

The output directory contains:

- `taylor_design.tsv`: one row per primary run, eligibility flags, family
  linkage, patient fold, and outcome/exposure availability flags;
- `patient_group_folds.tsv`: track-specific patient-to-fold assignments;
- `acquired_family_design.tsv`: longitudinal family completeness and readiness;
- `taylor_analysis_preparation.json`: counts and readiness/status summary.

Only metadata and availability indicators are written.  Raw AUC, GI50, and
cumulative-exposure values are intentionally absent from every output.

## Leakage controls

The fold assignment is deterministic and group-based on `patient_id`; random
OCM-level splitting is forbidden.  Feature selection, model fitting, and any
network/module scoring must be performed inside each training fold.  CORNETO
network inference remains label-blind and must not use outcome values to choose
priors, sources, targets, or modules.

