# Manuscript evidence ledger

This directory is the machine-readable evidence layer for `main.tex`.

- `claims.tsv` records each scientific claim, its estimand, a pre-specified
  falsification rule, current status, source, and wording boundary.
- `failures.tsv` preserves failed or superseded analyses that change scientific
  interpretation.
- `paper_ocm_evidence.tsv` distinguishes directly named OCM IDs from later
  crosswalks, count-only panels, non-OCM systems, and unknown exact sets.
- `roihu_result_snapshot.json` is a compact receipt-of-receipts. It contains
  selected values from scientific receipts on Roihu together with the original
  path, SHA256, JSON pointer, retrieval date, and claim limit. It is not a
  replacement for the original receipts.
- `numeric_ledger.tsv` gives each result number a source, locator, derivation,
  evidence class, unit, and claim limit.
- `study_ocm_registry.tsv` and `manuscript_evidence_snapshot.json` are generated
  by `scripts/build_manuscript_evidence.py` from frozen repository inputs.

Run:

```bash
python3 scripts/build_manuscript_evidence.py
```

Generated LaTeX is written to `tex/generated/`. Do not edit generated files by
hand. Scheduler state and status-document prose are never treated as scientific
receipts.
