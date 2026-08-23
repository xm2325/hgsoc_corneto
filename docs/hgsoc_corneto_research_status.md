# HGSOC CORNETO research status and dependency register

Last operational update: 2026-08-23 11:10 BST (13:10 EEST). This file is the project-level source of
truth for scientific scope, completed evidence, queued analyses, failed
attempts, dependencies, and claim limits. Slurm `COMPLETED` is never sufficient
on its own: a result is scientifically complete only when its output receipt
passes the corresponding content validator.

## Operational correction: 2026-08-20

The first checkpoint preparation jobs 727669, 727673, and 727677 failed before
writing context receipts because their three solver entry points did not export
the project WLS license file. Gurobi therefore selected a size-limited local
license and rejected the 8,400-row/26,192-column model. This was an execution
environment defect, not an OOM, infeasible model, or biological result. The
prepare, independent, and joint sbatches now export and readability-check the
same `GRB_LICENSE_FILE` used by the established solver jobs.

Corrected chains are 749575--749578 (E-MTAB-7223), 749579--749582
(E-MTAB-10801), and 749583--749586 (E-MTAB-14568). At the startup gate,
749575, 749579, and 749583 completed with exit 0 and wrote `prepared` context
receipts for 9, 13, and 27 samples. Logs confirmed WLS academic license
2849103. Their independent arrays remain guarded behind both
successful context preparation and the terminal state of active jobs 727583/727584, preserving the
solver-session safety margin. No scientific receipt exists yet, so none of
these jobs supports a biological claim.

## Timeout recovery update: 2026-08-23

Monolithic jobs 727583 (E-MTAB-11000) and 727584 (four-sample pooled smoke)
reached the 72-hour limit without writing scientific receipts. The pooled full
successor was therefore cancelled by dependency. This establishes operational
non-completion under the current monolithic formulation, not model
infeasibility and not a biological result; no further blind pooled retry was
submitted.

For each of E-MTAB-7223, E-MTAB-10801, and E-MTAB-14568, independent checkpoint
tasks 0 and 1 reached their original 24-hour limit without receipts, while
tasks 2 and 3 remained active. Slurm did not permit extending the running array
time limits. Idempotent 72-hour recovery arrays 805006, 805009, and 805012 were
therefore queued after the corresponding original arrays; they cover every
index but immediately skip any receipt that the original array completes. New
joint/assembly successors are 805007--805008, 805010--805011, and
805013--805014.

E-MTAB-11000 now has the missing checkpoint chain 805002--805005, with
independent concurrency one to preserve the WLS session margin. Fail-closed
availability, four-cohort comparison, and TPI1 successors are
805015--805018. No canonical cohort receipt exists yet.

## Central scientific question

> Do HGSOC OCMs contain metabolic and regulatory states that reproduce across
> cohorts, remain robust to repeated patients and modelling choices, are
> enriched relative to stroma/reference samples, and receive support from
> external mechanistic evidence?

The primary evidence unit is the 60 primary HGSOC tumour OCMs. They represent
52 patients and four studies, so analyses must distinguish OCM-level,
patient-balanced, cohort-stratified, and pooled results. The remaining public
RNA runs are reference/QC data and are not silently merged into the primary
cohort.

## Frozen sample universe

| Study | All RNA runs | Primary HGSOC tumour OCMs | Principal reference runs |
|---|---:|---:|---|
| E-MTAB-7223 | 36 | 9 | 16 stroma, 2 cell-line controls, non-HGSOC/ambiguous and other exclusions |
| E-MTAB-10801 | 36 | 13 | 17 stroma and 6 non-HGSOC |
| E-MTAB-11000 | 12 | 11 | 1 non-primary tumour run |
| E-MTAB-14568 | 33 | 27 | 6 non-HGSOC |
| **Total** | **117** | **60 OCM / 52 patients** | **57 reference/excluded runs** |

Reference categories overlap scientific roles and must be analysed separately:
stroma is the strongest formal secondary contrast; confirmed non-HGSOC is an
exploratory histology contrast; ambiguous samples are not controls; two
cell-line controls are descriptive QC only; later-passage and non-Tighe samples
support stability or response-blind replication rather than population-level
inference.

## Analysis hierarchy and claim limits

| Evidence tier | Analyses | Permitted interpretation |
|---|---|---|
| Primary | Four cohort-specific primary-HGSOC metabolic and regulatory models; cross-cohort recurrence | Reproducible model-predicted metabolic/regulatory state |
| Robustness | pooled-60, patient-balanced-52, lambda, PKN, NMF rank, alternative optima | Sensitivity to modelling and repeated patients |
| Reference | tumour vs stroma; confirmed HGSOC vs confirmed non-HGSOC; cell-line and passage QC | Enrichment/shared core, not absolute presence/absence |
| Mechanistic | Meeson/TPI1, WT vs delta-TPI1, FVA, order sensitivity | Model/known-mechanism concordance, not drug-response validation |
| Phenotype-linked | exact AUC, raw GI50, cumulative exposure, grouped CV | **Blocked until exact phenotype tables pass intake QC** |

All expression-derived fluxes are model-predicted feasible flux states, not
measured flux. Regulatory/NMF/metabolic agreement derived from the same RNA
profiles is internal consistency, not independent validation.

## Completed and auditable results

### RNA and pooled-input gates

- All four RNA studies were downloaded, quantified, and aggregated with strict
  run/receipt checks.
- Pooled primary expression contains 60 runs, 60 OCMs, 52 patients and 60,609
  genes; study counts are 9/13/11/27.
- Pooled metabolic input gate job **599942** completed and validated matrix,
  manifest, Human-GEM, and all SHA256 provenance. Receipt:
  `data/processed/corneto/pooled_primary_60/metabolic_input_gate.json`.

### NMF

- Primary cohort and pooled NMF jobs **591326-591328** completed.
- Pooled rank-3 clusters contain 23/20/17 samples; cophenetic correlation is
  0.939 and silhouette is 0.770. Rank 2 is more stable (0.972/0.899), so rank 3
  is a frozen comparability benchmark rather than a uniquely optimal rank.
- Pooled-vs-cohort rank-3 ARI is 0.723 (7223), 0.529 (10801), 0.377 (11000),
  and 0.569 (14568); mapped assignment agreement is 0.889/0.846/0.727/0.741.
- Study association with pooled state was not statistically clear
  (chi-square p=0.207, Cramer's V=0.265), but study-aware sensitivity remains
  required.
- Patient-balanced-52 rank-3 comparison completed: mapped agreement 0.942,
  ARI 0.842 and NMI 0.825. This supports robustness to repeated-patient
  weighting, but it is only one deterministic patient-balanced selection.

### Regulatory CORNETO

- True multi-condition normalized-lambda grid and retries completed; the final
  summary contains 45/45 pooled/cohort receipts across nine nominal lambdas.
- Patient-balanced regulatory analysis retained similar networks on the 52
  common patients: pooled/balanced union Jaccard 0.890 and mean per-sample
  Jaccard 0.899.
- Narrow-vs-richer PKN sensitivity completed. Pooled union Jaccard was 0.202
  and mean sample Jaccard 0.108; cohort union Jaccards were approximately
  0.108-0.143. Network conclusions are therefore materially PKN-sensitive and
  must be reported as stable cores plus uncertain alternatives.
- Regulatory longitudinal summary covered 60 runs and eight within-family
  transitions. It is response-blind; acquired-resistance interpretation still
  requires exact exposure and phenotype.
- Regulatory x NMF state integration completed for all four cohorts. It is a
  descriptive same-input integration, not independent validation or subtype
  discovery.
- E-MTAB-14568 regulatory alternative-optima ensemble and serial repairs
  completed. The solution ensemble, rather than a single sparse optimum, is the
  appropriate evidence object.

### Human-GEM and TPI1 model gate

- Public Meeson TPI1 table audit job **591424** completed; this is an audit of
  published/model values rather than an OCM-specific knockout result.
- Independent model-level TPI1 gate job **599873** completed and is valid:
  Human-GEM contains 13,096 reactions and 3,628 genes; TPI1 maps to
  `ENSG00000111669` and reaction `HMR_4391`; no solver was called. Receipt:
  `data/processed/corneto/tpi1_model_gate.json`.
- Obsolete pending job **591049** was cancelled before starting because it
  would have overwritten that same receipt after the original baselines. Its
  roles are now separated into 599873 (model only) and 599950 (strict
  receipt-dependent preflight).

## Terminal-job audit and retry lineage

This table records final evidence rather than counting every diagnostic Slurm
attempt as an independent experiment.

| Analysis family | Final evidence state | Failed/superseded attempts and disposition |
|---|---|---|
| RNA quantification and aggregation | Four aggregation receipts are `completed`: 36/36/12/33 runs, each with 60,609 genes and 227,462 transcripts | Earlier per-run download, OOM and aggregation failures were repaired; no RNA retry remains |
| Primary and pooled NMF | 591326-591328 completed; patient-balanced NMF/compare 592020 and 592094 completed | 592083 comparison failed and was replaced by 592094 |
| True multi-condition regulatory lambda grid | Final 591593 retry and 591595 summary completed; 45/45 receipts across pooled plus four cohorts and nine lambdas validate | Initial 591416 tasks hit Gurobi session limits; failed artifacts are preserved and were replaced only for affected labels |
| Regulatory alternative optima | Final 591572 summary validates 27 E-MTAB-14568 samples: 26 nonempty completed and one blocked zero-edge sample | Six initial session-cap errors were rerun serially by 591569 |
| Richer-PKN sensitivity | 592021/592023/grid jobs and comparison 592118 completed | No unresolved retry |
| Patient-balanced regulatory sensitivity | 592143 and 592149 completed on 52 common patients | No unresolved retry |
| Longitudinal regulatory summary | 592019 completed: 60 runs and eight response-blind within-family comparisons | Earlier prototype 588876 failed and was replaced by 588883/592019 lineage |
| Regulatory x NMF integration | Final v4 job 592053 completed with full 9/13/11/27 coverage; 576 BH-adjusted edge-state tests yielded no q<0.05 findings | 592018, 592040 and 592046 failed during interface/path corrections and were superseded by 592053 |
| Meeson public evidence | 591424 completed; toy joint/order/global-retention receipts validate algorithmic behavior | Cohort-specific order/ensemble jobs for 7223/10801/11000 were dependency-cancelled and have not yet been reattached to valid metabolic retries; 14568 tasks remain pending |
| Metabolic growth-fraction sensitivity | No valid receipt from 588286-588289 | All four reached the original 4 h limit; no retry is queued before primary baselines validate |

For the normalized regulatory grid, solver completion at high lambda is not a
positive network finding: pooled networks become empty at nominal lambda 0.05
and above, and most cohort networks also collapse at 0.05-0.1. This is evidence
of over-regularisation under the current scaling, not evidence of biological
absence. At lambda 0.001, pooled-vs-merged-cohort edge-union Jaccard is 0.746;
at lambda 0.01 it falls to 0.286. These remain response-blind technical results.

## Metabolic baseline: terminal failures and checkpointed recovery

Frozen scientific settings remain unchanged: Human-GEM v1.4.1, raw TPM transformed
with log1p, primary tumour only, candidate budget 25, growth fraction 0.9,
independent lambda 0.1, joint lambda 1.0, and explicit Gurobi with no fallback.

The original monolithic design solved every independent condition serially, then the
joint problem, and wrote the only receipt at process exit. Its terminal evidence is:

| Study | Original / retry evidence | Final receipt at 2026-08-19 audit |
|---|---|---|
| E-MTAB-7223 | 588250 TIMEOUT at 24 h; 600005 TIMEOUT at 72 h (step MaxRSS 52.6 GB) | missing |
| E-MTAB-10801 | 588251 OOM at 64G; 600004 TIMEOUT at 72 h/196G (step MaxRSS 60.5 GB) | missing |
| E-MTAB-11000 | 588252 TIMEOUT at 24 h; 600006 OOM after 43 h 30 min/128G | missing; targeted 196G retry **727583** is running |
| E-MTAB-14568 | 588253 OOM after 65 h 32 min/128G; 600007 TIMEOUT at 72 h/196G (step MaxRSS 57.4 GB) | missing |

Each retry log contains exactly the expected number of independent LP reads
(9/13/11/27) before termination. This supports an execution-stage diagnosis: completed
independent work was discarded while the process entered or attempted the joint stage.
It does not establish any biological result. Availability audit **599875** therefore
correctly reports status=incomplete, 0/4 valid receipts and
final_comparison_permitted=false; strict comparison **599876** failed closed on the
missing E-MTAB-7223 receipt. Dependent TPI1 jobs 599950/599951 were automatically
cancelled and produced no scientific output.

A checkpointed recovery was implemented without changing scientific parameters. A
frozen context records samples, candidate IDs, bounds, model/input SHA256 and solver;
each independent OCM writes an atomic receipt; the joint solve receives a separate
72-hour allocation; the canonical full_direct_b25.json is assembled only when all
condition receipts and the joint receipt are optimal and share the same context hash.
Targeted tests (checkpoint fail-closed assembly plus existing joint-FBA tests) pass 4/4.

| Study | Prepare | Independent array | Joint | Assemble |
|---|---:|---:|---:|---:|
| E-MTAB-7223 | 727669 | 727670 (0-8, %2) | 727671 (196G/72h) | 727672 |
| E-MTAB-10801 | 727673 | 727674 (0-12, %2) | 727675 (196G/72h) | 727676 |
| E-MTAB-14568 | 727677 | 727678 (0-26, %2) | 727679 (384G/72h) | 727680 |

The independent arrays are held until 727583 and pooled smoke 727584 terminate; their
combined maximum is six Gurobi tasks, below the project safety ceiling of eight active
solver workflows. Audit **727685** and strict comparison **727686** wait on the three
assemblies plus 727583 and remain fail-closed.

## Pooled metabolic joint analysis

The pooled analysis remains joint-only; it never repeats 60 independent MILPs. Both
previous four-condition one-per-study smoke attempts failed operationally without a
receipt: **600008** TIMEOUT at 12 h and **615508** TIMEOUT at 36 h. The latter reached
only 11.5 GB MaxRSS and showed no OOM, licence-session or Python exception, so one final
scientifically identical smoke **727584** was submitted with the maximum 72-hour wall
time and 128G. It is guarded by absence of the canonical smoke receipt.

If and only if 727584 succeeds, **727684** starts the pooled-60 joint solve (384G/72h).
Pooled-vs-cohort comparison **727689** requires both 727684 and strict four-cohort
comparison 727686. The earlier 600009/600010 and 615515/615517 jobs were cancelled by
failed afterok dependencies and never ran; they are provenance, not negative
scientific results.

## TPI1/FVA and Meeson dependency chain

The independent model-only TPI1 gate remains valid: Human-GEM contains 13,096
reactions and 3,628 genes; TPI1 maps to ENSG00000111669 and reaction HMR_4391.
It does not depend on cohort receipts and is not an OCM knockout result.

The scientific chain is now strict comparison **727686** -> receipt-dependent preflight
**727687** -> WT versus delta-TPI1 FBA/FVA **727688**. Any invalid cohort receipt cancels
this chain. Cohort-specific Meeson order/ensemble analyses also require the corresponding
canonical cohort receipt and will not be released from scheduler state alone.

## Reference comparisons still to implement/run

1. **Tumour vs stroma (formal secondary analysis).** Analyse E-MTAB-7223 and
   E-MTAB-10801 separately, use paired contrasts where identity permits, then
   meta-analyse prevalence/effect sizes. Report tumour-enriched,
   stroma-enriched, shared stable core, study-specific and unresolved features.
2. **HGSOC vs confirmed non-HGSOC (exploratory).** Confirm histology first;
   ambiguous samples remain a separate category.
3. **Tumour vs two cell-line controls (descriptive QC).** No population-level
   p-values or claims are justified.
4. **Passage/non-Tighe stability.** Use paired/case-level expression, NMF,
   regulatory and metabolic retention; do not attach drug-response labels.
5. For every network contrast, replace simple set subtraction with prevalence
   difference, patient/study-aware uncertainty, alternative-optima frequency,
   lambda/PKN sensitivity and multiplicity control.

## External response-blind validation

External results are auditable in
/scratch/project_2012997/xiaomei/hgsoc_corneto_external/regulatory_external_comparison_v2.json
(status=completed). The Taylor signature was frozen before external scoring at nominal
lambda 0.001: an edge required prevalence in at least 6/52 Taylor patients and 2/4
cohorts. Five signed edges passed; signature SHA256 is
e10fb081217ad15ece5119ed003dbc16f76cdf3699a1b78ee7bce41a6cf558e6.

True multi-condition fits were optimal for all evaluated conditions: 22/22 GSE277107
ovary/omentum bulk samples, 59/59 GSE189955 patient pseudobulk conditions, and 14/14
GSE208216 organoid models. At the prespecified 50% patient/model prevalence threshold,
none of the five Taylor edges was a consensus feature in any external group (ovary
n=11, omentum n=11, HGSOC epithelial-candidate n=12, fibroblast proxy n=12, normal-FT
n=6, PDO n=11, or FT organoid n=3). This is negative falsification evidence for the
strict external-consensus claim, not proof that the pathways are absent. Point
prevalences and bootstrap intervals remain descriptive; no drug-response or causal
claim is permitted.

GSE180661 source acquisition completed: the 32,555,276,423-byte HDF5 has locally
observed SHA256 edc5f5d7449478d7dfec1f575c89670bcd1bb041124ef2d53a4df0ecf7e29be6.
Because GEO supplied no upstream digest, the receipt correctly keeps
matrix_identity_frozen=false; pseudobulk remains blocked until this observed hash is
reviewed and frozen. DepMap likewise remains deliberately blocked until one explicit
quarterly release and its matching Model, Chronos and README files are supplied.

## Missing data and blocked claims

Exact OCM-linked AUC, raw GI50 units, cumulative treatment exposure and their
primary provenance are not yet available in a validated intake table. Until
they pass the frozen phenotype gate, the project cannot claim Taxol response,
intrinsic resistance, acquired resistance prediction, or response-model
accuracy. `chemo_naive_at_biopsy` is not an acceptable replacement for exact
cumulative exposure.

## Operational update rules

1. Preserve failed logs and receipts; retries write only when the canonical
   result is absent, and existing receipts must validate before a retry skips.
2. Use `afterany` for audit jobs so failures remain visible; use `afterok` plus
   internal receipt validation for scientific downstream jobs.
3. Record job ID, exact command, model/input hashes, sample contract, solver,
   lambda, candidate budget, objective and claim limit in every receipt.
4. Do not select lambda, NMF rank, PKN, or candidate budget based on the most
   attractive biological story. Primary settings and sensitivity grids are
   frozen before outcome analysis.
5. Update this register whenever a job becomes terminal or a validated receipt
   changes the evidence state.

## Recurring monitor

The Codex heartbeat **HGSOC CORNETO Roihu pipeline monitor** is attached to this
conversation and runs at a low-frequency approximately two-hour cadence. It has no explicit model or reasoning
override, so it follows the current conversation/default settings and does not
create a separate standalone monitoring conversation. This cadence is deliberate:

- Slurm dependencies already launch valid successors, so minute-scale polling
  would not accelerate the pipeline.
- Thirty minutes is short enough to catch smoke/startup, license-session,
  timeout, and receipt failures while 588250/588252 approach their 24-hour
  limits.
- It is long enough to avoid repeatedly querying multi-hour Gurobi jobs whose
  logs are sparse during branch-and-bound.

Each run checks scheduler state, resource use, log tails, and receipt JSON. It
may make only deterministic in-scope repairs: no duplicate submissions, no
cancellation of healthy jobs, no scientific-parameter changes, and no result
claim without a valid receipt. OOM/TIMEOUT retries preserve parameters and
increase only justified resources; Gurobi session-cap failures wait until fewer
than eight sessions are active. Meaningful state changes are written here and
pushed to the delivery branch; unchanged checks produce only a compact
checkpoint. After all cohort, pooled, comparison, and TPI1/FVA outputs are
terminal and audited, the monitor reports completion and should be paused.
