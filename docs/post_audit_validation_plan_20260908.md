# Post-audit validation: execution contract, 2026-09-08

## Question and authorized scope

Do expression-derived OCM states transfer to held-out patients/studies, and
do metabolic constraints actually use sample information under defensible
boundary assumptions? The user authorized the next validation stage and
biological interpretation of existing results. This does not authorize
releasing frozen b25 review holds or treating partial incumbents as biology.

The old nine running solver tasks are preserved. New work below is
licence-free (SciPy HiGHS LP and scikit-learn NMF), in the separate namespace
`data/processed/post_audit_validation_20260908/`. No old receipt is overwritten.

## Work packages and gates

| Package | Actual computation | Gate and claim limit |
| --- | --- | --- |
| Metabolic smoke | Four OCMs, one per study; null, frozen b25, real GPR scales 0.1/1/10 and two shuffled-gene controls; default vs exclusion of three energy uptakes | 15-minute job. All cases must be terminal LP outcomes with valid primal residuals for optimal solutions. Infeasible/zero-growth outcomes are preserved, not repaired by dropping constraints. |
| Metabolic information pilot | Eight OCMs, two per study; same policies with ten shuffles | Only after smoke receipt passes. At most one-hour job, no sparse MILP and no full-60 solve. |
| NMF smoke | Four patient-disjoint held-out-study folds, fixed ranks 2 and 3, train-only 500 MAD genes, two starts and two patient-balanced draws | 15-minute job; convergence and exact patient/run coverage checked. |
| NMF validation | Same folds, 6000 train-only MAD genes, 20 starts and 30 patient-balanced draws for each fixed rank | Only after smoke receipt passes; six-hour CPU job. Rank is not chosen on test results. Reconstruction versus training mean is a transferability diagnostic, not proof of subtype replication. |
| Existing biology | Hallmark annotation of actual rank-2/3 gene loadings, family-level expression distances, patient/study-matched tumour–stroma programme scores, actual selected-edge and imposed-output frequencies | Fifteen-minute independent job. Exploratory annotation, not clinical classification, measured metabolic flux or causal regulation. |

Full pilot and NMF continuation have both scheduler `afterok` and receipt
checks. Output directories must not exist before a run. Intermediate JSON is
saved after each case/fold; a failed/incomplete computation is not promoted.

## Metabolic assumptions that must remain visible

- Eight-sample selection is response-blind: lowest and median frozen cap
  counts within each study. It is a diagnostic panel, not the study denominator.
- Nested GPR AND uses minimum and OR uses sum. Missing constituents propagate
  unknown capacity. Zero RNA is not blindly converted to a knockout. Mixed
  AND/OR is evaluated, unlike the historical reproduction helper.
- Positive `log1p(TPM)` capacities and factors 0.1/1/10 are **uncalibrated
  sensitivity assumptions**, not measured enzyme capacity or flux units.
- Gene-value shuffles retain the marginal distribution and zeros among
  measured Human-GEM genes, disrupting the gene/expression association.
- `EX_atp[e]`, `EX_pep[e]`, `EX_pcreat[e]` uptake exclusions test dependence on
  those model permissions. **This is not an OCMI medium reconstruction.**
- Alternative scenarios maximize biomass without imposing the old unvalidated
  90% growth floor. They diagnose information/boundary dependence; they cannot
  be substituted into old canonical downstream gates.

The published OCMI recipe provides a real starting point: Nelson et al.,
Methods “Cell culture” and “Establishment of ex vivo models”
([source](https://www.nature.com/articles/s41467-020-14551-2)). The abbreviated
source record is `configs/ocmi_medium_evidence_20260908.json`. Exact basal
product variants, full ingredient/exchange mapping, per-study applicability
and uptake-rate assumptions remain unresolved; concentrations and gas
fractions are not uptake fluxes. No calibrated OCMI claim is released.

## Statistical and biological interpretation rules

- OCM74 spanning studies is excluded from training whenever any of its
  samples is held out. Training has one response-blind OCM per patient.
- Factor contribution labels normalize NMF basis scale. Repeated-draw ARI
  versus historical consensus is a multi-change sensitivity comparison,
  not a pure estimate of patient-selection variance or new validation cohort.
- NMF annotation uses the top 200 positive component-specific loadings after
  column-sum normalization; background is the mapped selected-gene universe,
  not the whole genome. BH covers all rank/state/Hallmark tests together.
- Hallmark version is fixed at 2025.1.Hs from the official Broad download;
  URL and checksum are in the output receipt. No GMT is staged in Git.
- Tumour–stroma uses within-sample mean percentile-rank gene-set scores,
  first averages samples within patient/study, then averages within patient.
  Family sign-flips and patient bootstrap intervals are exploratory and rely
  on their stated assumptions. These scores are not GSVA, GSEA or flux.
- Same-patient distances have matched study-pair references and equal-weight
  reference patient pairs. Shared reference pairs prevent naive independence;
  report descriptive effects without a spurious OCM-pair p-value.
- Regulatory target prominence must be shown alongside imposed input/output
  frequencies. A MYC-directed network may arise partly from requiring MYC as
  an output, not from discovering MYC independently.
- NMF/Hallmark/regulatory use overlapping RNA information. Agreement is
  internal consistency; tumour–stroma contrasts retain lineage/culture effects.

## Still pending after this execution stage

1. Complete and validate the OCMI reconstruction; then test numerical settings
   and fixed-indicator feasibility in a bounded sparse CORNETO/iMAT pilot.
2. Rebuild regulatory preprocessing separately within each training fold.
   Reusing the globally standardized existing bundle would leak held-out data;
   the current package does not claim a held-out regulatory refit.
3. Separate PKN breadth, source/output selection and depth one factor at a time.
   No new large Gurobi network grid is launched while the nine old sessions run.
4. Only after these gates: full metabolic analysis, joint/independent comparison,
   TPI1/FVA and external perturbation validation. Existing Meeson and Taylor
   findings are context, not a substitute for new OCM-specific validation.

## Execution receipts

Job identifiers and audited outcomes are recorded in the bilingual project
status register after submission. Submission alone is never success.

Submitted jobs: NMF smoke/validation 1121178/1121179; metabolic smoke/pilot
1121180/1121181; biological annotation 1121182. See
`evidence/post_audit_validation_submission_20260908.json` for observed evidence
and deployed source hashes. The metabolic smoke failed before LP execution at
the duplicate gene-ID guard. A local normalization repair preserves `_PAR_Y`,
but its deployment and the exact remote collision diagnosis are pending.
Remote reads were interrupted by approval-service usage limits. Do not infer
new numerical findings from completion logs or resubmit without checking successors.

Local verification after the normalization repair: all 12 tests in
`tests/test_post_audit_validation.py` passed; scoped Ruff, sbatch `bash -n`,
and `git diff --check` passed. These checks validate code behavior, not the
unread remote scientific outputs. No licence or downloaded GMT is included
in the delivery set; the unrelated uncontrolled LP JSON remains excluded.
