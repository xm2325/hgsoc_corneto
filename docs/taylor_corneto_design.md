# Taylor/CORNETO analysis design checkpoint

This document records the current feasibility boundary for the Taylor-group
paclitaxel question.  It is deliberately separate from the Meeson metabolic
benchmark and does not use any result from `csf3_slurm`.

## Scientific estimand

The primary question is whether label-blind CORNETO signalling modules add
patient-level information about paclitaxel response beyond the following
baselines:

| Model | Features available before response is used |
| --- | --- |
| M0 | ABCB1 and prespecified clinical/treatment covariates |
| M1 | gene expression |
| M2 | inferred transcription-factor activity |
| M3 | CORNETO modules |
| M4 | Barnes state, ABCB1, and CORNETO modules |

Paclitaxel AUC is the primary outcome.  `log10(GI50)` is secondary.  Response
values must not enter source selection, prior-network construction, CORNETO
optimization, or stability tuning.  Predictive evaluation is grouped by
`patient_id`; an OCM-level random split is invalid because several OCMs belong
to the same patient.

The three estimands remain distinct:

* **Cross-sectional:** between-OCM association using contemporaneous baseline
  tumour RNA and response.
* **Intrinsic:** the same association restricted to explicitly zero-exposure
  material.  `chemo_naive_at_biopsy` is retained as a proxy, not as proof of
  zero cumulative paclitaxel exposure.
* **Acquired:** within-patient changes in network features versus response
  changes, interpreted with intervening treatment exposure.  The claim limit
  is treatment-associated rewiring, not treatment-caused rewiring.

## Current cohort scaffold

The four public BioStudies accessions contain 117 paired-end RNA runs (82
tumour, 33 stroma, and 2 cell-line controls).  Applying the frozen HGSOC,
tumour-only, representative-library policy gives 60 eligible tumour runs from
60 OCMs and 52 patients:

| Accession | Eligible tumour runs |
| --- | ---: |
| E-MTAB-7223 | 9 |
| E-MTAB-10801 | 13 |
| E-MTAB-11000 | 11 |
| E-MTAB-14568 | 27 |
| **Total** | **60** |

Only 7 of these 60 rows are marked chemotherapy-naive at biopsy.  This is a
small sensitivity subset and remains a proxy until exact cumulative exposure
is available.  The manifest currently contains no exact paclitaxel AUC,
GI50, or cumulative-exposure values; therefore no response association has
been fitted or interpreted.

## Acquired-resistance feasibility

The current public RNA manifest contains multiple eligible tumour libraries
for seven patients:

`OCM66 (1/5)`, `OCM74 (1/3/5)`, `OCM110 (1/9)`, `OCM288 (4/7)`,
`OCM296 (3/5)`, `OCM327 (1/3)`, and `OCM333 (1/3)`.

The Taylor anchors `OCM231-1/5` and `OCM341-1/3` are not both represented in
the current four-accession RNA manifest (only `OCM231-1` and `OCM341-1` are
present).  They cannot support a within-patient RNA change until a matching
second library and phenotype are identified.  `OCM361a/b` remain explicitly
excluded from the longitudinal estimand by the analysis contract.

## Analysis gates

1. Finish Salmon quantification and one-reference aggregation for the three
   remaining accessions; retain all runs for QC and subset the primary cohort
   only after the manifest join.
2. Audit mapping, library-size, gene-detection, accession/batch, and duplicate
   passage fields before any signalling inference.
3. Run a small response-blind CORNETO signalling stability pilot.  The pilot
   tests source policy, prior-network choice, regularization, and module-level
   recurrence; it is not a search for a single “resistance edge.”
4. Add the exact Taylor phenotype table when it is available, then fit the
   prespecified M0–M4 comparisons with patient-grouped validation.  If M3 does
   not improve held-out performance, that negative result is the valid answer.
5. Treat acquired analyses as exploratory until paired response and exposure
   values exist for the same patient/timepoint pairs.

The Meeson metabolic/CORNETO benchmark is a separate track and should not be
used to fill the missing Taylor phenotype fields.
