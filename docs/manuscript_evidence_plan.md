# HGSOC CORNETO manuscript evidence and verification plan

## Objective

Build `main.tex` as both a falsifiable scientific manuscript and a navigable
project wiki. Every quantitative statement must be traceable to one of four
evidence classes:

1. **Observed source fact** — copied from a primary paper, BioStudies/ENA
   record, frozen supplementary table, or scientific receipt, with an exact
   locator and checksum where available.
2. **Deterministic derivation** — recomputed from a versioned input by a tested
   script; the derivation and denominator are explicit.
3. **Model output** — reported only after a receipt-level scientific gate, with
   solver, parameters, status, and claim limits.
4. **Configured or planned analysis** — clearly labelled as such and never
   written as an observed result.

Scheduler state, log existence, status-document prose, and a non-empty output
file are not scientific evidence by themselves.

## Pyramid structure

The manuscript is organised from decision-relevant conclusions to supporting
detail:

1. **Answer first:** what is supported, weakened, falsified, blocked, or still
   pending for the central HGSOC OCM proposition.
2. **Primary evidence:** frozen cohort construction, NMF concordance,
   patient-weighting sensitivity, regulatory robustness, optimisation
   ambiguity, and the current metabolic/phenotype boundary.
3. **Methods and falsification:** estimands, nulls, perturbations, receipt
   gates, and observations that would overturn each claim.
4. **Wiki appendices:** study/run/OCM registry, paper-to-OCM registry,
   claim ledger, failed/superseded analyses, public resources, and data-access
   boundaries.

## Build tasks

### A. Freeze the evidence schema

- `evidence/claims.tsv`: claim, estimand, falsification rule, status,
  evidence locator, permitted wording, prohibited wording.
- `evidence/failures.tsv`: scientifically informative failed or superseded
  analyses, root cause, correction, and residual claim limit.
- `evidence/paper_ocm_evidence.tsv`: paper-level reported counts, exact OCM IDs
  when recoverable, experimental role, source locator, and uncertainty.
- `evidence/roihu_result_snapshot.json`: selected receipt values, original
  Roihu paths, SHA256, JSON pointers, retrieval date, and claim limits.

### B. Generate all derived cohort numbers

`scripts/build_manuscript_evidence.py` must regenerate:

- the 117-run study/OCM/patient registry;
- Table 1 counts and exclusion categories;
- the 60-OCM/52-patient repeated-patient arithmetic;
- FASTQ counts and exact byte totals;
- LaTeX macros and long-table rows;
- a machine-readable evidence snapshot with hashes of every local source.

Generated files are committed because arXiv does not run project-specific
Python during compilation. They must never be hand-edited.

### C. Rewrite and shorten the header

- Keep the full scientific title on page 1.
- Set `\shorttitle` to `Cross-cohort HGSOC OCM states` so the arXiv running
  header does not contain the full title.
- Move detailed resource material to appendices while preserving direct
  navigation from the main text.

### D. Replace non-falsifiable prose

- Split the compound central proposition into independently testable claims.
- Remove the inference that low edge overlap proves stable modules.
- Describe rank 3 as a frozen comparison anchor, not the optimal biological
  rank; report rank-2 stability as a counter-result.
- Describe pooled-versus-cohort NMF as same-data concordance, not independent
  external replication.
- Describe one lexicographic 52-patient analysis as one deterministic
  sensitivity, not general patient-level robustness.
- Keep metabolic, TPI1 OCM dependency, tumour-reference enrichment, and drug
  response claims blocked until their specified evidence gates pass.

## Forward verification plan: source to manuscript

1. Rebuild all generated evidence files from frozen inputs.
2. Confirm expected invariants:
   - 117 unique runs;
   - 60 primary HGSOC tumour OCMs;
   - 52 unique patients;
   - study split 9/13/11/27;
   - exclusions 33 stroma, 2 controls, 17 non-HGSOC or ambiguous tumour, four
     unmatched-to-Tighe tumour, and one non-representative duplicate;
   - 234 FASTQ objects and 477,762,645,114 compressed bytes.
3. Verify every remote model number against the recorded receipt SHA and JSON
   pointer.
4. Confirm each paper-to-OCM row is either direct, a later authoritative
   crosswalk, count-only, non-OCM, or unknown; never silently infer an ID.

## Reverse verification plan: manuscript to source

1. Extract every numeral and percentage from `main.tex` and generated TeX.
2. Classify each as citation metadata, source fact, deterministic derivation,
   configured parameter, model output, or date/version.
3. Require an evidence-ledger entry, source citation, equation, or explicit
   derivation for every scientific number.
4. Search for prohibited wording from the claim ledger and fail the check if it
   appears outside a quoted warning.
5. Check that `pending`, `blocked`, and `falsified` claims are absent from the
   positive-result summary.

## Structural and release checks

- Unit-test cohort classification, repeated-patient arithmetic, LaTeX escaping,
  invalid claim statuses, and deterministic output.
- Run the generator twice and require a clean second diff.
- Run `git diff --check` and Python tests.
- Check balanced LaTeX braces/environments, bibliography keys, duplicate labels,
  missing `\input` files, and forbidden Unicode where relevant.
- Compile the manuscript when a TeX engine is available; otherwise record the
  absence and run structural checks without claiming a successful PDF build.
- Stage only the manuscript, evidence ledgers, generator, generated TeX,
  tests, and this plan. Never stage raw data, solver licences, unrelated status
  changes, or scratch outputs.
- Push without force to the existing delivery branch; if it diverges, create a
  narrowly named manuscript branch.

## Acceptance criteria

The revision is acceptable only if:

1. Table 1 and the complete study-ID appendix are generated from the frozen
   manifest rather than typed by hand.
2. Every modelled value has a receipt path, checksum, and JSON/TSV locator.
3. Every substantive claim states what observation would falsify it.
4. Unknown exact OCM panels remain explicitly unknown.
5. Failed or superseded analyses are preserved in the appendix without being
   presented as scientific successes.
6. The paper distinguishes description, same-data robustness, independent
   validation, and causal evidence.
