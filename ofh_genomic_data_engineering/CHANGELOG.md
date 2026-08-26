# Changelog

## 0.4.0 - 2026-08-26

- Add a provider delivery manifest and a fail-closed delivery gate before genomic processing.
- Validate expected source SHA-256, declared reference genome, VCF sample count and deterministic sample-roster hash.
- Add deterministic delivery fingerprints and bind them into provenance and release parameters.
- Add idempotency semantics: an exact registered duplicate is `NOOP`; the same delivery ID with changed content is `REJECT`.
- Add six delivery-validator tests covering acceptance, checksum mismatch, wrong genome build, sample-roster mismatch, exact duplicate and delivery-ID collision.
- Expand deterministic resume validation to the 12-process graph and include delivery validation in tracked release-facing hashes.
- Synchronise Python and Nextflow manifest versions to 0.4.0.

## 0.3.1 - 2026-08-26

- Add an explicit deterministic rerun contract to GitHub Actions.
- Snapshot release-facing output SHA-256 values, rerun the same workflow with Nextflow `-resume`, and require every one of the 11 processes to be served from cache.
- Fail validation if any tracked BGEN, Parquet, schema, query, provenance or release output changes after the resumed run.
- Emit `reproducibility_validation.json` as workflow evidence and document the scope of the reproducibility claim.

## 0.3.0 - 2026-08-26

- Convert PLINK text outputs into typed Parquet instead of all-string tables.
- Store genomic positions and counts as `int64`; QC fractions, probabilities and PCA scores as floating point.
- Use ZSTD compression and emit a machine-readable Arrow schema manifest.
- Add DuckDB region-query validation cross-checked against pandas.
- Make the release stage depend on a successful query-layer contract and include schema/query evidence in provenance.

## 0.2.0 - 2026-08-26

- Add a source inventory derived from the VCF instead of provider naming.
- Add a release contract that checks schemas, row counts, sample preservation, variant-key uniqueness, QC value ranges and product hashes.
- Add deterministic release IDs and a Nextflow release gate that exits non-zero on failed data contracts.
- Add negative tests for duplicate sample IDs and tampered product hashes.
- Document the boundary between FAIR-compatible release metadata and GA4GH API implementation.

## 0.1.1 - 2026-08-25

- Handle PLINK2 `.pvar` VCF-style metadata before the tabular header when producing Parquet data products.
- Add regression tests for metadata preambles and missing PLINK headers.
- Move Nextflow parameter defaults into `nextflow.config` to remove undefined-parameter warnings.

## 0.1.0 - 2026-08-25

- Add the initial real-data genomic engineering workflow: pinned 1000 Genomes VCF, bcftools normalisation, PLINK2 import/QC/PCA, BGEN export, Parquet outputs, provenance, Docker and GitHub Actions validation.
