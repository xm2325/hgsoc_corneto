# HGSOC CORNETO research status and dependency register

Last execution update: 2026-09-09 (BST); latest bounded LP audit: 2026-09-09. This file is the project-level source of
truth for scientific scope, completed evidence, queued analyses, failed
attempts, dependencies, and claim limits. Slurm `COMPLETED` is never sufficient
on its own: a result is scientifically complete only when its output receipt
passes the corresponding content validator.

## 2026-09-10: exact-matrix witness and presolve diagnostic

User requested continued numerical repair and then withdrew the14:53 scheduling
request after quota reset. The attempted one-time scheduling update was rejected;
no new wake-up was created. Existing two-hour monitoring is unchanged.

New diagnostic **1240622** retains the same restricted model, growth floor and
gap1e-4, but disables HiGHS presolve (120s limit). Before solving, the runner
evaluates the saved all-indicators-on witness against the exact assembled
constraint matrices and variable bounds and writes residuals to receipt.
Any nonvalidated outcome now exits nonzero while preserving the receipt.
Updated binary regression test and Ruff pass. Output is a fresh
`data/processed/revised_binary_support_20260910_nopresolve` directory;
no old evidence or hold is overwritten. Outcome pending; presolve is a tested
hypothesis, not yet an established cause.

## 2026-09-10: pFBA support gate passed; restricted binary test submitted

Follow-up:1237008 is scheduler COMPLETED but receipt is no_incumbent,
HiGHS infeasible (status2),767 restricted reactions, no objective/bound/gap.
Receipt SHA256 `b43d0aba47b630a946b027cac435cbddb16bc2bd130e7eb71a46c87a199a1e20`.
This fails the binary validation gate despite the corresponding continuous
fixed-support LP passing. Diagnose exact matrix/bounds against the known feasible
flux before attributing this to presolve or changing scientific constraints.
No global solve is authorized by this result. Medium gate remains unresolved.

1232423 completed in21s. All eight fixed-support rechecks passed; four closed-
uptake negative controls had zero biomass. Input and panel hashes matched;
maximum pFBA mass residual1.308e-12, bound violation0. Receipt SHA256:
`355769d1af982bdebd16ec4244a25e52ff790e180576aee0b51478c22ff573d3`.
Support counts in receipt order:816,767,827,934,809,887,792,869. These are
thresholded pFBA supports, not proven minimum-cardinality networks.

New job **1237008** tests binary links on the first sample's90%-growth pFBA
support only: minimize indicator count with per-reaction L*y<=v<=U*y,
120s HiGHS MILP limit, requested relative gap1e-4, saved incumbent/bound/gap,
full-primal/link/integrality checks and an independent fixed-support LP.
It is topology-restricted and not a full CORNETO result. No Gurobi session,
10min/32G allocation;15 local tests passed. Output:
`data/processed/revised_binary_support_20260910/receipt.json`.
Medium calibration and unrestricted binary model remain unfinished gates.

## 2026-09-10: user-authorized revised-model smoke submitted

New job **1232423**,4CPU/32G/15min, tests full GPR scale1 caps with pFBA,
50%/90% of each sample's own maximum growth, four fixed response-blind OCMs,
and four all-boundary-uptake-closed negative controls. No Gurobi sessions.
Flux-selected support must pass a second LP with every unselected reaction zeroed.
Fourteen local tests and scoped Ruff passed. Fresh output:
`data/processed/revised_gpr_pfba_20260910/receipt.json`.

This is a continuous parsimonious-flux diagnostic, not minimum-cardinality
CORNETO or a calibrated OCMI model. Medium reconstruction remains a release gate;
old b25 holds remain intact. No full-cohort successor was submitted. Monitor
1232423 as newly authorized work; audit its cases/negative controls and hashes
before any further action. See [revision contract](revised_model_validation_20260910.md).

## 2026-09-10 08:52 UTC: legacy terminal accounting and gate evidence

Targeted sacct confirms joint/assembly 834324–834331 and TPI1/FVA
834334–834335 CANCELLED without execution. Comparison 834333 FAILED because
the 7223 full_direct_b25.json was absent. Although 834332 is scheduler COMPLETED,
its log and receipt explicitly say incomplete (log: ready=0); this is not a
successful four-cohort analysis. Availability receipt SHA256:
`92ae4f477bad8ada484dfb6a54811ebd8a42fc89d396044827f8474ec25bf8e0`.

Example rapid failure 863034_0 is an intentional scientific_review_hold rejection,
not OOM. Example 70-hour failure 937737_3 ends user_limit/partial_incumbent;
receipt `003_ERR6389080_job1077075_task3.json` SHA256 is
`fb2a7b56aad8251d7635ba5a238d7dca72b348a706bbb3fa59d9a2f2bb715dfa`.
These examples do not classify every array failure. Remaining individual terminal
artifacts need review; do not infer cause solely from exit code or duration.
No retries, cancellations or hold releases were performed. Bounded pilot remains
completed; its numerical validity does not repair the old sparse-network chain.

## 2026-09-09 08:20 UTC: expression-shuffle contrast checked

Pilot receipt SHA remains unchanged. At GPR scale1 under legacy default bounds,
all eight real-expression maxima exceed all ten within-sample gene-value shuffles.
Real maxima span 6.2809–10.2168; shuffle maxima across all 80 cases span
0.00850–4.24632. This supports gene/value identity mattering under this particular
GPR-cap construction, not just the marginal expression distribution. It is not
HGSOC specificity, measured growth, patient discrimination or a validated medium.
Only ten shuffles per sample were run: a corrected one-sided Monte Carlo tail
estimate is 1/11 when none exceeds the real value, not p=0. Shared seeds and
gene/model structure preclude treating the 80 cases as independent patients.

Targeted squeue returned no entries for the documented legacy chain members.
That is queue absence only; terminal accounting/artifact review is still pending.
No retries or new research jobs were submitted; scientific holds remain in place.

## 2026-09-09 03:20 UTC: bounded pilot completed and numerical receipt audited

1163554 completed in 1m58s. Receipt: 240 unique LP cases across eight samples,
all optimal and all numerical_pass. Eight input and four context SHA256 checks,
model hash and both code hashes match. Maximum mass-balance residual is
6.884e-11, maximum bound violation 5.685e-13. These are audited continuous-LP
diagnostics, not canonical sparse CORNETO networks.

In both tested boundary settings, expression-null and frozen-b25 growth maxima
remain identical across eight samples at 187.35362997658078. Real GPR caps yield
sample-dependent maxima: scale0.1 0.6281–1.0217; scale1 6.2809–10.2168;
scale10 30.4407–35.0413. This demonstrates sensitivity to these assumed caps,
not measured metabolic capacity. Excluding the three energy uptakes does not
change these maxima (apart from rounding); thus their presence alone does not
explain the growth optimum. A full medium audit remains necessary. Shuffle
contrasts and legacy terminal-artifact review remain pending; no holds released
or further full-scope jobs submitted. See `evidence/metabolic_pilot_audit_20260909.json`.

## 2026-09-08 22:08 UTC: serialization failure repaired; two-hour monitoring

Post-deployment check: 1163513 COMPLETED in 22 seconds; its receipt reports
completed, 56 LP cases across four samples and zero unresolved solver cases.
Pilot 1163554 is PENDING (Priority), with smoke dependency cleared. Detailed
policy contrasts and source/hash re-audit remain pending; no new mechanism claim.

Targeted sacct/log audit: smoke 1138529 FAILED after 14 seconds; pilot 1138550
CANCELLED without running. The PAR_Y guard was passed (3,627/3,628 model genes
covered for the first sample). One LP result was saved, then NumPy's boolean
from the bound-violation comparison could not be JSON encoded. The stale receipt
still says running: scheduler/log evidence takes precedence for execution state,
and this partial diagnostic is not a completed scientific result.

Explicit conversion of numerical_pass to Python bool is deployed. A regression
test simulates a small nonzero bound residual and verifies JSON serialization;
all 13 tests and scoped Ruff checks passed. Replacement smoke **1163513** and
pilot **1163554** use fresh `_r3` directories; pilot requires afterok:1163513
and the matching smoke receipt. No scientific parameters or legacy holds changed.
Deployed runner SHA256: `e3dd9399b9e10efb753e2daa9bb9d2b2ee1e1eb51f22b10f841eb5559644ffb0`.

The existing same-conversation heartbeat now checks every two hours, without
model override or a separate task. It reads this register for direct successors,
stays quiet on unchanged/non-actionable state, and reports meaningful failures,
audited results or required action. Old upload prohibitions were removed because
authorized GitHub delivery through 98feba1 is verified.

## 2026-09-08 14:31 UTC onward: connection restored and receipts checked

SSH succeeded on retry; the prior approval-service quota block is historical,
not current. `sacct` confirms 1121178/1121179/1121182 COMPLETED, 1121180 FAILED,
1121181 CANCELLED. Full NMF receipt is completed, all eight recorded input hashes
match, all eight NPZ artifacts exist, and all eight fold patient overlaps are empty.
Biology receipt is completed and its ten top-level source hashes match. This does
not constitute independent validation of the regulatory subsection's nested sources.

Full NMF patient-mean held-out reconstruction improvement vs training mean:
7223 rank2/rank3 12.06%/13.85%; 10801 10.62%/14.79%; 11000 9.34%/14.26%;
14568 8.92%/12.62%. Thirty patient-balanced draws per rank yield ARI median
0.9225 (rank2, range 0.5255–1.0000) and 0.8708 (rank3, 0.8421–0.9468).
These are transferability/sensitivity diagnostics, not validated clinical subtypes;
rank3's extra capacity also improves reconstruction. Six of seven repeated families
have lower within-patient than matched between-patient expression distance;
OCM327 is the exception. Paired tumour–stroma tests cover 17 patients, not 60;
29/50 Hallmark tests have BH q<0.05. Lineage/culture and shared RNA remain limitations.

Remote input inspection confirmed 45 first-dot collisions, each exemplified by
ordinary/PAR_Y pairs; version-only normalization leaves zero collisions in the
checked common gene list. The fix is deployed. Replacement smoke **1138529**
and gated pilot **1138550** use fresh `_r2` outputs; pilot requires afterok:1138529
and the new smoke receipt. Old failed evidence is preserved. No b25 hold released.
See `evidence/post_audit_validation_audit_20260908.json` for hashes and provenance.

## 2026-09-08: bounded validation submitted; receipt audit interrupted (historical)

The user authorized the next validation stage. Five license-free jobs were
submitted using `hpc/roihu/post_audit_validation.sbatch`; none releases the
b25 application review holds or cancels healthy running work.

| Job | Purpose | Last directly observed evidence, not a current live claim |
|---|---|---|
| 1121178 | NMF smoke, 15 min | Log reports eight folds and completed receipt; receipt content not yet audited |
| 1121179 | Patient-grouped NMF validation, 6 h | PENDING (Priority), dependency cleared |
| 1121180 | Metabolic information smoke, 15 min | Failed before LP: version stripping creates duplicate gene identifiers |
| 1121181 | Metabolic information pilot, 1 h | Submitted afterok:1121180; terminal state not yet verified |
| 1121182 | Biological annotation, 15 min | Log reports completed receipt; results not yet read or audited |

The next remote read was rejected by the approval service's usage limit before
SSH executed. This is not evidence of a certificate failure. No remote repair,
replacement submission, or new GitHub delivery is claimed after that rejection.
The existing monitor schedule and prompt have not yet been updated for these IDs.

Local repair replaces lossy first-dot gene-ID truncation with Ensembl-version-only
normalization preserving `_PAR_Y`; genuine normalized collisions still fail closed.
The exact conflicting remote IDs remain unverified, so this is a tested code repair,
not proof of the particular input cause or a deployed fix. Do not overwrite the
failed smoke receipt. Before one replacement smoke, audit 1121181 and any successor,
use a fresh output namespace, and ensure its pilot reads that replacement receipt.

See [execution plan](post_audit_validation_plan_20260908.md) and
[biological interpretation and falsification limits](existing_biology_interpretation_20260908.md).
New NMF/annotation output must not yet be quoted as a validated biological result.

## 2026-09-07: authorized delivery completed; no Slurm hold applied

Following explicit user authorization, all 21 reviewed files, including the four
audit JSON files, were pushed to `xm2325/hgsoc_corneto`, branch
`agent/metabolic-retry-pooled-status`, at commit `3907c4f` (fast-forward from
`e7f03b7`). The four JSON files were copied to the original Roihu project's
`evidence/` directory; all four remote SHA-256 hashes match the local files.
This supersedes the historical delivery-blocked notes below. No raw data,
licences, secrets, or the untracked uncontrolled LP audit were included.
The user requested clarification of Slurm hold; no scheduler-level hold was
applied. Existing application startup review holds remain unchanged. This
delivery is not a new scientific result or authorization for new solver work.

## 2026-09-07 08:31-08:36 BST: read-only continuation, no new scientific result

The one-time 04:00 BST resumption was delayed: its heartbeat timestamp was
05:16 BST and the first executed check was 08:31 BST. The cause was not
established. Using the existing same-conversation automation and
[official OpenAI documentation](https://learn.chatgpt.com/zh-Hans/docs/automations),
the prior **Thursday 09:00 Europe/London** recurrence
was restored and verified, with no model override or new task.

Targeted `squeue` still shows the same nine RUNNING tasks (about 38-41 h).
All nine log tails are current. Live gaps at the 08:35 check are:

| Array tasks, in order | Actual numeric JobIds, in order | Rounded live gaps |
|---|---|---|
| 937737_3 / 4 / 6 | 1077075 / 1077163 / 1077866 | 3.11% / 3.62% / 3.59% |
| 948765_3 / 4 / 5 | 1077699 / 1077700 / 1077798 | 4.01% / 3.94% / 3.83% |
| 834322_11 / 12 / 13 | 1076518 / 1077046 / 1077074 | 3.07% / 3.94% / 3.14% |

`sstat` on these actual `.0` steps reports MaxRSS samples of 23.44-33.81 GiB
and nonzero CPU/I/O counters; this is not a guarantee against later OOM.
`sacct` remains unavailable because its database connection was refused, so
no new terminal accounting conclusion is asserted.

The four context hashes and all 20 saved attempt hashes match the prior audit;
there are no new attempts, canonical independent files, or joint receipts.
All four application review holds remain `hold`. The 11000 `afterany`
serialization and downstream dependencies remain unchanged. No solver was
submitted, modified, cancelled, or placed on Slurm hold.

Read-only GitHub verification still returns `e7f03b7d630b33d48dc16fcbaa3c5f30c07067aa`.
Prior upload/Slurm-hold denials are not new authorization: neither blocked
transfer nor hold was retried. This checkpoint is committed **locally only**;
GitHub, OneDrive and remote copies do not yet include this checkpoint.

Latest operational follow-up: **2026-09-06 23:02 BST**. Targeted `squeue`
returned the same nine RUNNING tasks; all nine solver-log tails were current.
Live gaps were 3.14/3.65/3.62% for 937737_3/4/6,
4.05/3.97/3.87% for 948765_3/4/5, and 3.10/3.97/3.17% for
834322_11/12/13. These rounded values are not final receipts or completion
estimates. [Saved log-tail evidence](../evidence/scientific_audit_running_logs_20260906.json).
An attempted scheduler-level hold of pending solver jobs was rejected by the
safety reviewer before execution. **No Slurm hold was applied**; it needs
explicit user approval. The previously deployed application startup hold
remains in place. No running job was targeted or changed.

**Delivery follow-up, 23:09 BST:** implementation/audit commit `fa7efa7` exists
locally in `/private/tmp/hgsoc-scientific-audit-20260906`. The GitHub push to
`xm2325/hgsoc_corneto`, branch `agent/metabolic-retry-pooled-status`, was
rejected before execution by the safety reviewer; it is **not uploaded**.
The four new audit-evidence JSON transfers back to Roihu were separately
rejected. Do not bypass either rejection or retry without explicit payload
and destination approval. Reviewed solver code and the bilingual documents
were deployed; JSON evidence remains available locally. OneDrive is not
claimed synchronized. The 04:00 continuation must read this isolated working
copy or the updated remote documents, not assume GitHub already has this audit.

## Latest decision: b25 scientific review hold, 2026-09-06

The research question remains worthwhile, but continuing identical long b25
solves is not currently justified as OCM-specific biological analysis. See the
[full scientific audit and repair report](hgsoc_corneto_scientific_audit_20260906.md),
[60-context/receipt audit](../evidence/scientific_validity_audit_20260906.json),
and [controlled fixed-indicator LP audit](../evidence/fixed_indicator_lp_controlled_audit_20260906.json).

- All 60 contexts share growth optimum 187.35362997658078; six have no
  expression-derived caps. All 20 saved independent attempts (17 distinct
  OCM/runs) report zero flux on their applied expression-capped reactions.
  Each saved flux satisfies 50-60 of the 60 context-bound sets; 15 satisfy all
  60. This is weak specificity evidence, not proof of identical feasible sets.
- All 20 saved indicator selections are **infeasible** in a continuous LP at
  the original growth floor. All 20 positive controls retaining the union of
  indicator-selected and reported nonzero-flux reactions are **feasible**.
  Consequently these selections cannot be used as growth-supporting networks
  for knockout/FVA. Sparse-summary mass balance alone did not detect this.
- Canonical independent files remain **0/9, 0/13, 0/11, 0/27**. Nineteen
  70-hour attempts ended at gaps 2.9679%-4.1531%; the additional 600-second
  smoke has gap 5.8292%. These are optimization diagnostics, not biology.
- The 14 historical snapshot source hashes all still match. The 45 regulatory
  grid receipts all report optimal solver status, although 27 have empty edge
  unions. None has a square edge-by-condition shape affected by the repaired
  export bug. This does not certify every other historical regulatory output.

**Repairs are deployed and backed up on Roihu.** Canonical acceptance now
requires solver/gap, full-primal and artifact-hash validation. Joint startup
checks every independent receipt in the cohort, not only a successful subset
repair array. Regulatory partial incumbents no longer masquerade as completed
fits. The changes passed 29 targeted tests and remote Python syntax checks;
[deployment evidence](../evidence/scientific_audit_deployment_20260906.json)
records code/context hashes and four successful pre-solver hold checks.

Four `checkpoint_b25/scientific_review_hold.json` markers now block **new**
independent/joint starts before solver import. This is an application startup
hold, **not a confirmed Slurm administrative hold**. No healthy RUNNING process
was cancelled or changed. Such deliberately blocked starts must not trigger
automatic retries. No frozen context, MIPGap or medium parameter was altered.

The last successful targeted queue snapshot showed `937737_3/4/6`,
`948765_3/4/5` and `834322_11/12/13` running. Later Slurm controller calls timed
out, so that nine-task snapshot is not a real-time guarantee. SSH and file
audits remained available. E-MTAB-11000 serialization was corrected and
verified as `afterany:1083050_*`, `afterany:1083051_*`, `afterany:863034_*`;
unrelated cohort success is not a scientific prerequisite for 11000. The
joint/assembly/comparison/TPI1 fail-closed gates remain intact.

Next priority is a versioned medium/GPR/expression-policy and numerical pilot
with null/shuffled-input controls, not another identical 70-hour retry.
Patient-grouped held-out NMF/regulatory validation can be developed in parallel.
The dated operational history below is provenance, not authority to release
the review hold or recreate old jobs.

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

**Interpretation hold, documented 2026-08-28:** the saved b25 partial solutions
expose sparse sample-specific expression constraints, uncalibrated exchange
bounds and flux/indicator numerical disagreement. Eventual solver convergence
alone will not resolve these issues. The frozen b25 runs remain optimization
benchmarks; patient-specific metabolic claims require a separate input and
numerical-quality audit. No scientific parameter was changed by this monitor.

### Within-patient coverage, verified 2026-08-27

Filtering `evidence/study_ocm_registry.tsv` on `primary_cohort_eligible=true`
and grouping by `patient_id` gives seven repeated patients with 15 primary
OCMs: six two-OCM families and one three-OCM family. Together with 45 singleton
patients this reproduces 60 OCMs and 52 patients. One baseline-to-each-other
contrast gives eight comparisons, not eight independent patients.

| Patient | Primary OCM IDs | Already-submitted independent array tasks |
|---|---|---|
| OCM66 | OCM66-1; OCM66-5 | 834320_5; 834320_4 |
| OCM74 | OCM74-1; OCM74-3; OCM74-5 | 834320_7; 834320_6; 834322_26 |
| OCM110 | OCM110-1; OCM110-9 | 834321_2; 834321_1 |
| OCM288 | OCM288-4; OCM288-7 | 834322_4; 834322_5 |
| OCM296 | OCM296-3; OCM296-5 | 834322_8; 834322_9 |
| OCM327 | OCM327-1; OCM327-3 | 834322_15; 834322_16 |
| OCM333 | OCM333-1; OCM333-3 | 834322_18; 834322_19 |

The 14568 indices also have the existing 863034 repair coverage. These are
independent sample solves, not a newly submitted patient-level joint model.
Jobs 834324-834327 remain cohort-level joint models. Regulatory job 592019's
`regulatory_longitudinal_joint_l0p001.json` was re-read: it records eight
response-blind comparisons in these seven families, not treatment causality.
Its baseline labels do not independently establish collection chronology.

Patient-level joint sensitivity was discussed but not launched. A test of
within-patient similarity needs independent solutions and matched
between-patient comparisons; joint union regularization already favours
reaction reuse and cannot itself establish a patient effect. OCM74 spans
7223 and 14568, so cohort-specific candidate selection is an additional
confound. Chronology, treatment and same-biopsy/spatial relationships must be
verified from source metadata rather than inferred from suffix ordering.

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
- Narrow-vs-richer graph-policy sensitivity completed. Pooled union Jaccard was 0.202
  and mean sample Jaccard 0.108; cohort union Jaccards were approximately
  0.108-0.143. Network conclusions are therefore materially PKN-sensitive and
  must be reported as stable cores plus uncertain alternatives. Inputs,
  outputs and depth changed together, so this is not a single-factor test of
  PKN breadth.
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

## Metabolic baseline: historical active/failed/queued snapshots

The dated snapshots below retain retry provenance. The 2026-09-06 scientific
review hold at the top of this document supersedes their continuation rules.

Primary settings are frozen: Human-GEM v1.4.1, raw TPM transformed with
`log1p`, primary tumour only, candidate budget 25, growth fraction 0.9,
independent lambda 0.1, joint lambda 1.0, and explicit Gurobi without fallback.

All monolithic cohort attempts are terminal without a valid final receipt.
Jobs 588250/588252 and 600004/600005/600007 reached their wall-time limits;
588251 failed after a genuine 64G OOM; 588253 and 600006 failed without a
canonical scientific receipt. The later 11000 retry 727583 also reached 72 h.
These attempts established that a final-only JSON write is not an adequate
checkpointing strategy for this MILP.

The replacement design freezes each cohort's expression-derived bounds and
objectives in `checkpoint_b25/context.json`, solves one independent OCM per
array task, solves the joint cohort only after every independent receipt is
canonical, and assembles `full_direct_b25.json` last. At this update:

| Study | Required OCM receipts | Previous r4/legacy outcome | Cancelled queued 64G tasks | Current recovery array |
|---|---:|---|---|---:|
| E-MTAB-7223 | 9 | 805860_5: partial incumbent at 70 h; others OOM | none remained pending | **834320**, 128G, tasks 0-2 running; 3-8 queued |
| E-MTAB-10801 | 13 | 805861_1: partial incumbent at 70 h; other started tasks OOM | 805861 tasks 6-12 | **834321**, 128G, tasks 0-2 running; 3-12 queued |
| E-MTAB-11000 | 11 | 805003_0: Slurm TIMEOUT at 72 h, no receipt | 805862 tasks 0-10 | **834323**, 128G, `0-10%3`, waiting for the other three r5 arrays |
| E-MTAB-14568 | 27 | r4 started tasks OOM | 805863 task 8 and tasks 12-26 | **834322**, 128G, tasks 1/3/5 running, 6-26 queued; r5 tasks 0/2 failed OOM and 4 failed SIGBUS; **863034**, 256G, pending repair |

At the 00:24 EEST live check, the six 24-hour legacy tasks had run for
13 h 13 min and the 72-hour 11000 task for 11 h 44 min. Solver-step CPU
efficiency was approximately 94-96%, disk counters continued to increase, and
peak RSS was 15.2-23.8 GiB against 64G requests. This supports active solving
rather than OOM or an idle hang, but it does not establish feasibility,
optimality or scientific correctness. These pre-instrumentation processes have
no live Gurobi progress log and have written no canonical or partial receipt,
`.sol` or `.mst`; their incumbent and gap are therefore unavailable. The six
24-hour tasks have about 10 h 47 min remaining and may still reproduce the
earlier timeout pattern. Instrumented jobs provide live solver logs and write
end-of-limit telemetry/partial receipts when the solver returns normally;
an OS kill or failed filesystem write can still prevent persistence.

The six healthy running tasks from 749576/749580/749584 and the healthy task
805003_0 were not cancelled. Earlier tasks 0-1 in each of
749576/749580/749584 had already timed out at 24 h without independent
receipts. Each instrumented recovery task first validates and skips any
matching canonical receipt, so completed legacy work is not recomputed.

After the instrumented arrays started, repeated r4 elements across 7223, 10801
and 14568 were killed by the Slurm memory cgroup, commonly within 1-9 minutes,
at the original 64G request. The failed Python steps did not get an opportunity
to write canonical or partial scientific receipts. At 12:09 EEST, seven
healthy 64G tasks remained running and were deliberately retained. Every still
pending 64G task was cancelled, including the entire not-yet-started 805862
array. The independent-job default was raised to 128G and recovery arrays
**834320-834323** were accepted with throttle three. They validate and skip
matching canonical receipts, so only missing/OOM indices are recomputed.

The r5 dependencies use the explicit Slurm array wildcard (`jobid_*`), rather
than only the array master ID, so the recovery arrays cannot start while a
retained r4 task is still running. The 11000 array additionally waits for the
three other recovery arrays and legacy job 805003, keeping the planned solver
load within the ten-session Gurobi licence ceiling. This is an operational
recovery only: there is still no newly validated biological result.

At the 03:42 EEST check, r4 tasks 805863_10 and 805863_11 had also ended in
`OUT_OF_MEMORY` after about 8 h 6 min; the existing 834322 recovery already
covers both indices, so no additional retry was submitted. Five instrumented
r4 tasks remained healthy and actively solving: 805860_5, 805861_0,
805861_1, 805861_5 and 805863_0. Their live relative gaps were 3.30-3.91%.
Legacy 805003_0 remained active at 43 h 30 min but has no live telemetry.
No cohort had a canonical independent receipt or joint receipt, so these
incumbents remain optimization diagnostics only. All r5 arrays correctly
remained blocked on their explicit dependencies.

At the 17:11 EEST check, 805861_0 and 805863_0 had subsequently ended in
`OUT_OF_MEMORY` after approximately 27 h and 30 h, respectively. Their existing
128G recovery arrays cover both indices; no duplicate retry was added. With
the last 14568 r4 task terminal, 834322_0 started normally and produced an
initial feasible incumbent with a 17.7% relative gap after about 12 seconds;
the remaining 14568 tasks were pending for scheduler priority, not dependency
or licence failure. The still-running r4 tasks 805860_5, 805861_1 and 805861_5
had live gaps of 3.26-3.87%. Legacy 805003_0 remained active at 52 h 32 min.
Canonical independent and joint receipt counts remained zero in every cohort,
so no biological interpretation was released.

Task 834322_0 then failed with a second cgroup `OUT_OF_MEMORY` event after
5 min 42 s despite its 128G request. It had reached an incumbent of
-134.95360, bound -143.05700 and 6.00% relative gap immediately before the
kill, but wrote no canonical or partial receipt; those values are diagnostic
only. Repair array **863034** was therefore submitted for all 27 indices with
256G, throttle three, the same frozen context and unchanged solver parameters.
It waits for all 834322 tasks, validates and skips matching canonical receipts,
and recomputes only missing/noncanonical indices. The 14568 joint job 834326
now depends on `afterok:863034_*`; this removes the permanently unsatisfiable
dependency on the failed r5 array while preserving fail-closed progression.

At the 09:50 EEST check, 805861_5 had also ended in `OUT_OF_MEMORY` after
31 h 5 min, and 834322_2 had ended in `OUT_OF_MEMORY` after 4 h. Existing
recovery arrays 834321 and 863034 already cover the respective indices, so no
new successor was required. Active solver tasks 805860_5, 805861_1,
834322_1, 834322_3 and 834322_4 remained healthy with live gaps of
3.39-3.83%; legacy 805003_0 remained active at 69 h 11 min without live
telemetry. The remaining r5 14568 tasks were limited by the array throttle,
and r6 remained correctly blocked until r5 terminates. Canonical independent
and joint receipt counts remained zero across all four cohorts.

### 2026-08-27: audited full-duration partial outputs

Two instrumented r4 tasks returned at their 252000-second solver limit and
persisted their results before Slurm termination. Both receipts explicitly
report `partial_incumbent`, `scientific_success=false`, Gurobi `TIME_LIMIT`
and matching frozen-context SHA256 values. Their Slurm `FAILED 2:0` is the
intentional fail-closed exit for partial output, not a new OOM.

| Array task / RNA run | Incumbent objective | Best bound | Relative gap | Receipt suffix under the cohort's `checkpoint_b25/instrumented_attempts/` |
|---|---:|---:|---:|---|
| 805860_5 / ERR2808261 | -136.0536082251 | -141.1906456601 | 3.775745% | `005_ERR2808261_job833065_task5.json` |
| 805861_1 / ERR6389069 | -136.2536081782 | -141.1262291334 | 3.576141% | `001_ERR6389069_job832834_task1.json` |

Each receipt's referenced `.sol`, `.mst` and `.gurobi.log` exists and is
nonempty; `summary_error` is null. This validates result persistence, not
optimality or biological conclusions. Both gaps remain above the unchanged
`MIPGap=0.0001` (0.01%) threshold. The deployed runner and reviewed code hashes
match: the runner writes `.mst` files but does not reload them for a retry.
Thus recovery currently reuses canonical completed samples, not the partial
incumbent or a saved branch-and-bound tree. No warm-start reuse is claimed.

Existing 834320/834321 arrays started automatically and already cover both
partial tasks; no additional retry was submitted. Legacy 805003_0 reached
72 h without a receipt and is covered by 834323. At this snapshot nine 128G
solver tasks were active: three each in 834320, 834321 and 834322, with live
gaps 3.19-4.65%. Actual-step `sstat` reported approximately 5.2-29.6 GiB RSS
and nonzero CPU/I/O counters; these container accounting samples are not a
guarantee against OOM. All four cohorts still have zero canonical independent
and joint receipts. The audit, comparison and TPI1/FVA gates stay closed.
The same-conversation monitor now targets the r5/r6 chain explicitly, retaining
its existing schedule and no model override.

### 2026-08-27 partial-solution scientific audit, registered 2026-08-28

This read-only audit used the two named 70-hour attempt receipts above, their
`.sol` and `.mst` files, frozen `context.json` inputs, the deployed objective
and indicator code, and Human-GEM v1.4.1. The model SHA256 remains
`57d1b137f0c90d83a3e4f9a8225d74d37523594e6ee99f622b160a014d9f7050`.
Context SHA256 values are
`7ea9d2268ec7647bfa0f47f8215913442cafc2b6f76faefa13a40c402b7fcb1b`
(7223) and
`1475eee9397af6644fcfaa6500594fbe1441b7d4d4acfbf714a66bac91f0792c`
(10801). These are model diagnostics, not validated biological findings.

| Saved-solution quantity | OCM66-1 / ERR2808261 | OCM110-9 / ERR6389069 |
|---|---:|---:|
| Expression-derived reaction bounds actually applied, excluding biomass | 15 | 1 |
| Selected indicators, threshold >=0.5 | 513 | 511 |
| Nonzero fluxes, absolute value >1e-7 | 544 | 542 |
| Biomass flux | 187.3536299766 | 187.3536299766 |
| Nonzero fluxes with indicator below 0.5 | 31 | 32 |

- The cohort candidate budget is 25, not 25 effective expression bounds in
  every OCM. `_candidate_sets` in `scripts/run_corneto_14568_pilot.py` retains
  only candidates with `proposed_upper > 0`; `_reaction_bounds` skips missing
  candidates. Zero-expression candidates therefore do not automatically close
  reactions. Missing and zero expression require distinct treatment before a
  revised biological analysis; neither justifies an unreviewed zero bound.
- All expression-capped reactions are zero at the saved flux reporting
  threshold. Each saved flux vector also satisfies the other OCM's expression
  and biomass bounds at tolerance 1e-7. The vectors pass the SBML mass-balance
  check with maximum absolute residual below 3.34e-9 and have no model-bound
  violation above 1e-6. This establishes reciprocal flux feasibility, not
  equality of the complete feasible sets or uniquely patient-assigned networks.
- Both solutions import ATP, phosphoenolpyruvate and phosphocreatine through
  `EX_atp[e]`, `EX_pep[e]` and `EX_pcreat[e]`, each at -1000. The unchanged
  model allows these exchanges; the supplied context does not calibrate them
  to measured uptake or actual culture medium. Biomass and energy-pathway
  conclusions cannot be interpreted as measured OCM physiology.
- The 31/32 flux-indicator discrepancies have flux magnitude approximately
  0.00386-0.00906 and binary values approximately 3.86e-6-9.06e-6 in `.mst`.
  With the deployed `lb*y <= v <= ub*y` constraints and bounds up to 1000,
  these are consistent with integer-tolerance-amplified trickle flow. Passing
  the flux mass balance does not validate a network after rounding indicators.
  Tightened justified bounds and a fixed-indicator feasibility audit are
  needed before classifying active reactions. See the
  [Gurobi IntegralityFocus explanation](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#integralityfocus).
- The deployed single-sample objective is `-biomass + 0.1*sum(indicators)`.
  Recomputing it from `.sol` matches the reported incumbent objectives; the
  approximate 0.2 difference is explained by sparsity, not different biomass.
  The indicator sets share 402 reactions out of a 622-reaction union
  (Jaccard 0.6463), but this is not an HGSOC-specific core. PPP reactions
  `G6PDH2c`/`PGLc` and TPI1 reaction `HMR_4391` are used in these solutions;
  their presence is not evidence of essentiality, enrichment or drug response.

Input/media revision, expression-null controls, integrality checks,
alternative-solution/FVA analysis and patient-balanced contrasts remain
unperformed follow-up gates for these interpretations. Lowering the MIP gap
alone cannot substitute for them. No media, expression policy, lambda,
integrality tolerance or running job was changed during this audit.

### 2026-08-28 11:15-11:17 BST: live recovery and covered SIGBUS failure

Targeted `squeue`/`sacct`, current solver-log tails and `sstat` on the actual
numeric solver-step IDs establish the following snapshot. Gaps below are
rounded live-log values, not final receipts or estimates of completion time.

| Array task | Actual numeric JobId | Elapsed at inspection | Live gap |
|---|---:|---|---:|
| 834320_0 | 890339 | 27 h 54 min | 3.06% |
| 834320_1 | 890340 | 27 h 54 min | 3.72% |
| 834320_2 | 890341 | 27 h 54 min | 3.37% |
| 834321_0 | 890279 | 28 h 3 min | 3.47% |
| 834321_1 | 890280 | 28 h 3 min | 3.65% |
| 834321_2 | 890281 | 28 h 3 min | 4.19% |
| 834322_1 | 862558 | 68 h 1 min | 3.15% |
| 834322_3 | 864522 | 66 h 35 min | 3.42% |
| 834322_5 | 900491 | 21 h 19 min | 3.27% |

All nine logs were recently updated. Actual-step RSS samples were about
18.0-45.6 GiB with nonzero CPU/I/O counters; as before, container accounting
does not prove that a 128G allocation cannot OOM. The two oldest 14568 solves
are approaching their internal 70 h limits, not guaranteed convergence.

Newly audited failure **834322_4**, actual JobId **866509**, run
**ERR13907041 / OCM288-4**, ended on 2026-08-27 15:57:28 EEST after
42 h 38 min 29 s. Accounting reports `FAILED 7:0`; the terminal output
`logs/met-inst-ind-14568-r5-834322_4.out` explicitly reports a Singularity
wrapper **Bus error** and srun exit 135. This confirms SIGBUS; neither OOM nor
a Gurobi session-cap error is established. The underlying cause remains
unresolved. The last solver-log values were incumbent -136.25361, bound
-141.20775 and gap 3.64%. No canonical/attempt JSON, `.sol` or `.mst` was
written for that attempt, so only progress-log evidence survives. The existing
write-on-return instrumentation does not protect against abrupt process loss.

Array **863034**, already requesting 256G and throttle three, covers this
missing/noncanonical index as well as the earlier OOM indices. Its dependency
is still `afterany:834322_*`, and joint **834326** still requires
`afterok:863034_*`. The larger allocation is not a demonstrated fix for
SIGBUS. No duplicate or additional retry was submitted and no healthy job was
cancelled. All other joint, assembly, comparison and TPI1/FVA gates retain
their existing dependencies. The 11000 array **834323** remains pending on
the three r5 independent arrays using explicit wildcards.

Canonical independent receipts remain **0/9, 0/13, 0/11 and 0/27**; all four
joint and `full_direct_b25.json` files are absent. The two earlier 70-hour
partial receipts still match their contexts and retain their artifacts; there
is no new partial receipt from the r5 arrays at this snapshot. No biological
interpretation or new patient-level/pooled research job was released.

### 2026-08-29 12:05-12:08 EEST: node-local failures and covered 10801 repair

The targeted audit found nine active solvers and still no canonical independent
receipt (**0/9, 0/13, 0/11 and 0/27**). The active 7223 tasks 834320_0-2 had
run for about 50 h 43 min with live gaps 3.00%, 3.59% and 3.28%. Active 10801
tasks 834321_7/10/11 had run for about 4-6 h with gaps 3.86%, 4.00% and 4.24%.
Active 14568 tasks 834322_5/6/7 had run for about 19-44 h with gaps 3.21%,
3.82% and 4.38%. All nine logs were current, and actual solver-step RSS samples
were approximately 5.7-40.1 GiB with nonzero CPU and I/O counters.

Three long 10801 tasks (834321_0-2; actual JobIds 890279-890281) failed on
node `rc5140` with Singularity-wrapper SIGBUS after about 44-47 h. Their last
logged gaps were 3.43%, 3.61% and 4.14%; none wrote an attempt or canonical
receipt. Six immediately following tasks (834321_3-6 and 834321_8-9) also
landed on `rc5140` and failed in 7-10 seconds because `srun` could not execute
`/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python`, reporting
`No such file or directory`. This common-node pattern is operational evidence
for node/filesystem/runtime failure, not OOM, solver infeasibility or a Gurobi
session-cap error.

One fail-closed 10801 repair array, **937737** (`0-12%3`, 128G), was therefore
submitted with unchanged frozen context and solver parameters. It waits on
`afterany:834321_*`, excludes `rc5140`, validates and skips any future
canonical receipt, and recomputes only noncanonical indices. Joint job 834325
was rewired to `afterok:937737_*`. `rc5140` was also excluded from the still
pending elements and downstream solver jobs in the documented recovery chain;
already-running healthy tasks were not changed or cancelled.

For 14568, 834322_1 and 834322_3 both reached the 252000-second solver limit
and then ended with wrapper SIGBUS. Task 834322_1 wrote no attempt receipt.
Task 834322_3 atomically wrote and was audited as `partial_incumbent`,
`scientific_success=false`, matching context SHA256, Gurobi `TIME_LIMIT`,
objective -136.6536067185, best bound -141.3243116334 and relative gap
3.417916%; its referenced `.sol`, `.mst` and Gurobi log exist and are nonempty.
This is optimization evidence only. Existing 256G repair array 863034 covers
both noncanonical indices, so no duplicate 14568 retry was submitted. All
joint, assembly, comparison and TPI1/FVA scientific gates remain closed.

### 2026-08-30 12:00-12:03 EEST: 7223 partial receipts and repair serialization

The first three 7223 r5 tasks reached the unchanged 252000-second Gurobi limit.
Tasks 834320_0 and 834320_1 atomically wrote audited `partial_incumbent`
receipts before their Singularity wrappers ended with SIGBUS. Both receipts
have `scientific_success=false`, the matching 7223 context SHA256, Gurobi
`TIME_LIMIT`, requested `MIPGap=0.0001`, and nonempty `.sol`, `.mst` and solver
log artifacts. Their objective/bound/gap values are respectively
-137.0536041377/-141.1212280161/2.967907% and
-136.2536156425/-141.1159266378/3.568574%. Task 834320_2 also reached the
solver limit, with last logged objective -136.85361, bound -141.29013 and gap
3.24%, but SIGBUS occurred before an attempt receipt was written. None is a
canonical scientific result.

Because 834320 had no existing successor for its noncanonical indices, one
fail-closed 7223 repair array, **948765** (`0-8%3`, 128G), was submitted with
the same context and scientific/solver parameters. It waits on
`afterany:834320_*`, excludes the known-bad node `rc5140`, validates and skips
canonical receipts, and recomputes only noncanonical indices. Joint job 834324
now requires `afterok:948765_*`.

The 11000 independent array 834323 was serialized behind all three repair
arrays (`948765_*`, `937737_*` and `863034_*`). Thus at most three throttle-3
repair arrays, rather than those arrays plus 11000, can request Gurobi sessions
at once. This preserves the operational ceiling without changing scientific
parameters. Separately, 10801 task 834321_10 was killed by its 128G Slurm
memory cgroup after 16 h 55 min; the already-submitted 937737 repair covers it,
so no additional retry was added.

Nine solver tasks were active at this check. Live gaps were 4.25/4.18/4.24%
for 7223 tasks 3-5, 3.59/3.59/3.72% for 10801 tasks 7/11/12, and
3.18/3.58/4.28% for 14568 tasks 5-7. All logs were current; actual-step RSS
samples were approximately 6.0-53.9 GiB with nonzero CPU/I/O. Canonical counts
remain **0/9, 0/13, 0/11 and 0/27**. No joint, assembly, comparison or
TPI1/FVA scientific gate was released.

### 2026-09-06 13:31-13:35 EEST: partial-only progress and Slurm control outage

No canonical independent receipt exists in any cohort: counts remain
**0/9, 0/13, 0/11 and 0/27**; joint, assembly, comparison and TPI1/FVA outputs
remain absent. Since the previous checkpoint, 14 additional attempt receipts
were content-audited: five for 7223 (indices 3, 5, 6, 7 and 8), four for 10801
(r5 index 11 and r6 indices 0-2), and five for 14568 (indices 5, 7, 8, 9 and
10). Every receipt matches its cohort context, reports Gurobi `TIME_LIMIT`,
`scientific_success=false`, requested `MIPGap=0.0001`, has null summary error,
and references nonempty `.sol`, `.mst` and solver-log artifacts. Their gaps
span 3.18-4.15%. They are optimization evidence only and do not release any
scientific dependency.

Four r5 tasks without a receipt have explicit OOM evidence in their terminal
logs: 7223 index 4 after 99,206 solver seconds; 10801 indices 7 and 12 after
191,700 and 174,591 seconds; and 14568 index 6 after 251,025 seconds. Their
existing repair arrays cover the missing/noncanonical indices; no duplicate
was submitted.

At this inspection, nine solver logs were still advancing: 7223 r6 indices
3-5 (gaps 4.10/4.03/3.94%), 10801 r6 indices 3/4/6
(3.65/3.77/3.66%), and 14568 r5 indices 11-13
(3.31/4.03/3.21%). However, `squeue` and `scontrol` returned no job records,
and `sacct` reported that its persistent Slurm database connection was refused.
Four repair elements (7223 r6 indices 0-2 and 10801 r6 index 5) have terminal
logs stating that `srun` could not confirm their allocations because the Slurm
socket timed out and the JobIds were expired/invalid. This is scheduler/control
infrastructure failure, not solver or biological evidence.

Because scheduler state and allocations could not be authoritatively audited,
no new retry was submitted and no dependency was changed. The failed repair
indices will require one consolidated, fail-closed successor only after Slurm
control/accounting recovers and the currently running arrays can be audited.
The existing 11000 serialization and all downstream fail-closed gates were
left untouched.

### 2026-09-06 13:41-13:44 EEST: Slurm controller recovery and targeted repair

`squeue` recovered and showed the same nine solver tasks still running, while
`sacct` remained unavailable because the accounting database connection was
refused. No healthy solver was cancelled. To cover only the four confirmed
allocation-timeout failures, two fail-closed successors were accepted:

- **1083050**, E-MTAB-7223 r7, array `0-2%3`, 128G, excluding `rc5140`, after
  `afterany:948765_*`;
- **1083051**, E-MTAB-10801 r7, task `5%1`, 128G, excluding `rc5140`, after
  `afterany:937737_*`.

Both reuse the frozen cohort context and unchanged instrumented scientific and
solver parameters. They validate and skip a matching canonical receipt, so
they cannot overwrite valid work. Joint jobs **834324** and **834325** now wait
for `afterok:1083050_*` and `afterok:1083051_*`, respectively. The serialized
E-MTAB-11000 array **834323** at that checkpoint waited for both new repairs plus
`afterok:863034_*` (**subsequently corrected to afterany for all three parents
in the 2026-09-06 scientific audit above**). The dependencies were verified through `squeue`; the new
arrays remain pending behind their active parents and therefore add no current
Gurobi sessions. Canonical counts remain unchanged and no scientific result
was released.

The r5 instrumented independent tasks request 128G, eight CPUs and a 72 h Slurm limit,
but pass an internal Gurobi `TimeLimit=252000` seconds (70 h), `MIPGap=1e-4`,
eight threads and seed 0. This leaves time to atomically write an attempt
receipt before Slurm termination. On an incumbent it records objective, best
bound, absolute/relative gap, status, solution count, node/work/iteration
counts and model dimensions, and writes `.sol`, `.mst` and a solver log. A
time-limited incumbent is `partial_incumbent`, exits non-zero to block
`afterok`, and is explicitly **not** a biological result or a canonical cohort
receipt.

Smoke job **805824** validated the mechanism on E-MTAB-7223 OCM ERR2808250.
After exactly 600.0 s it had 10 feasible solutions, incumbent objective
-134.7536165, best bound -142.6087313 and relative gap 5.8292% after 14,986
nodes. It atomically wrote a `partial_incumbent` receipt plus `.sol` and `.mst`
artifacts, then intentionally exited 2. This is successful observability and
fail-closed control, not scientific completion. As of this update the four
cohorts still have 0/9, 0/13, 0/11 and 0/27 canonical independent receipts.

The replacement instrumented joint jobs are **834324-834327**; assembly jobs
are **834328-834331**. Joint memory is 196G for 7223/10801/11000 and 384G for
14568, with the same 70 h internal solver limit. Job 834326 is explicitly
gated by the 256G repair array 863034. Audit **805872** and strict
comparison **805873** were superseded and cancelled. New audit **834332** and
strict comparison **834333** wait on all four r5 assembly jobs with `afterany`
and will fail closed on missing or partial receipts. TPI1 preflight **834334**
and full WT-versus-delta-TPI1 FBA/FVA **834335** require `afterok` from the
strict comparison chain. The superseded r4 downstream jobs 805864-805875 were
cancelled only after this complete replacement chain had been accepted by
Slurm.

## Pooled metabolic joint analysis

This analysis is intentionally joint-only: the four cohort jobs provide
independent/cohort evidence, while repeating 60 independent MILPs inside one
72-hour pooled job would waste completed work and erase progress on timeout.

The pooled input gate **599942** remains valid, but all three final-only pooled
smokes failed operationally: 600008 timed out at 12 h, 615508 at 36 h and
727584 at 72 h. None wrote a scientific receipt; their full-60 and comparison
successors were cancelled without running. Repeating the same opaque pooled
job is therefore stopped.

The next pooled attempt must first build a frozen pooled checkpoint context and
use the same instrumented joint runner: a bounded two-to-four-condition smoke
must expose incumbent/bound/gap and write a partial receipt before a 60-condition
job is eligible. No new pooled solver job was submitted in this update because
the cohort checkpoint contexts cannot be silently substituted for a pooled
context. Pooled-vs-cohort biological interpretation remains blocked until both
the four canonical cohort receipts and an instrumented pooled receipt validate.

## TPI1/FVA and Meeson dependency chain

```mermaid
flowchart LR
    M["599873 model-only TPI1 gate: valid"]
    B["834333 four strict metabolic receipts"] --> P["834334 receipt-dependent TPI1/FVA preflight"]
    P --> V["834335 WT vs delta-TPI1 FBA/FVA"]
```

Existing Meeson order-sensitivity and ensemble jobs are attached to individual
original cohort jobs. If an original times out but its conditional retry
succeeds, the corresponding cancelled downstream task must be resubmitted
against the valid retry receipt. No result should be inferred from dependency
state alone.

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
6. Every long MILP must set an internal solver time limit shorter than Slurm,
   persist incumbent/bound/gap and solver artifacts, and distinguish
   `partial_incumbent` from canonical `completed` output. Partial receipts may
   guide recovery but may not release scientific downstream jobs.
7. The 2026-09-06 scientific review hold supersedes old retry rules: no
   automatic identical b25 long solves or retries of deliberate startup holds.
   A reviewed versioned input/numerical pilot must precede reconsideration.

## Recurring monitor

The Codex heartbeat **HGSOC CORNETO Roihu pipeline monitor** is attached to this
conversation without an explicit model/reasoning override or a separate task.
On 2026-09-06 the user-requested one-time resumption was configured for 16:32
Europe/London. The actual continuation signal arrived at 17:56:09 BST; the
reason for that delay was not established. Work continued in this conversation.
The prior actual recurrence, **Thursday 09:00 Europe/London**, was then restored
and verified in the automation configuration. The old text claiming a
30-minute recurrence was incorrect and has been removed.

**Historical scheduling update, 2026-09-06:** after the user reported that the
five-hour usage window had reset, the same heartbeat was rescheduled for a
one-time continuation at **2026-09-07 04:00 Europe/London (03:00 UTC)**.
This configuration was verified; execution at that future time is not yet
established at that checkpoint. The prompt resumes only unfinished work and restores the prior
Thursday 09:00 recurrence after that continuation. No model override or new
task was created.

The 2026-09-07 continuation above has now run, and the verified current
recurrence is Thursday 09:00 Europe/London. The one-time 04:00 rule is no
longer active. Delivery was subsequently explicitly authorized and completed
as recorded above. Administrative Slurm holds still await informed approval;
the monitor must not infer that approval from the upload authorization.

The revised prompt reads the scientific audit first, preserves healthy RUNNING
jobs, honours the four application startup review holds, and never retries
deliberately blocked starts or identical 70-hour b25 attempts automatically.
Only documented jobs/successors are inspected, without broad or high-frequency
polling. Controller outages do not establish terminal status. New frozen
model/media policies and biological full solves require a reviewed pilot and
authorization; existing context parameters and strict downstream gates remain
unchanged. Meaningful results/failures/actions update both documents and the
GitHub delivery branch; unchanged or non-actionable checks remain quiet. When
the preserved running work and authorized audits are terminal, report final
reviewed status and that this monitor can pause, without creating new scope.
