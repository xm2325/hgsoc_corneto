# Changelog

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
