# Taylor phenotype intake gate

This gate is intentionally phenotype-free until an exact, source-attributed
Taylor table is supplied. It does not download, OCR, infer, impute, transform,
merge, or analyse phenotype values. Missing input produces a blocked receipt
and exit code 2.

The external TSV is long-form, with exactly one row per primary OCM and
endpoint. Required columns are:

```text
canonical_ocm_id patient_id drug endpoint value unit endpoint_definition is_exact source_file source_record_id
```

Supported endpoints are `paclitaxel_auc` (primary), `paclitaxel_gi50`
(secondary; raw nM, with `log10` reserved for the later analysis stage), and
`cumulative_paclitaxel_exposure`. The source's exact AUC convention must be
recorded in `endpoint_definition`; the placeholder unit token
`source_reported_auc` avoids asserting that the published AUC is dimensionless.
Exposure is currently configured as `mg` to match the frozen manifest field.
If the exact source table reports another basis such as `mg/m2`, the schema
must be changed explicitly before intake. The validator never converts units.

Run only after the exact table is available:

```bash
PYTHONPATH=src python scripts/validate_taylor_phenotype_intake.py
```

The gate rejects duplicate OCM/endpoint records, proxy or imputed values,
unknown OCMs, patient-map conflicts, missing provenance, non-finite values,
invalid value domains, and unit mismatches. It evaluates completeness against
the frozen 60-OCM primary cohort. Its JSON receipt contains counts, missing OCM
IDs, readiness flags, and an input checksum, but no phenotype values.
Invalid input also writes an `invalid_phenotype_intake` receipt and exits 2, so
an intake failure cannot be mistaken for analysis readiness.

Association is permitted only when all 60 primary AUC rows pass. Secondary
analysis additionally requires all 60 GI50 rows. Exposure-stratified analysis
additionally requires all 60 exact exposure rows. A blocked or invalid gate is
not an authorisation to run partial association.
