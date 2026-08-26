# Data Release Contract

## Purpose

The pipeline treats a genomic delivery and a research data release as two separate contracts. Processing success alone is not enough to publish a release.

## 1. Source contract

After VCF normalisation, `source_inventory.py` records:

- the number of sample IDs read directly with `bcftools`;
- a SHA-256 digest of the ordered sample-ID list;
- the number of indexed variants;
- a SHA-256 digest of the normalised VCF.

Duplicate or zero sample IDs and zero indexed variants are rejected before release validation.

This avoids relying on filenames or provider-side labels. The current public fixture is a useful example: its filename includes `N100`, while the file itself contains 90 sample IDs.

## 2. Curated release contract

`validate_release.py` checks the BGEN, Parquet, summary and provenance products as one release.

The gate validates file presence, required columns, row-count consistency, unique sample IDs, unique normalised variant keys, source-to-release sample preservation, variant-count non-inflation, bounded QC metrics, cross-table sample consistency and SHA-256 integrity.

The gate emits `release_validation.json`. A passing release contains a deterministic `release_id` computed from source identity, QC parameters, summary metadata and product hashes.

## 3. Rejection path

The command exits non-zero when any check fails. The Nextflow `RELEASE_GATE` process therefore fails, and GitHub Actions does not mark the build as releasable.

Unit tests include explicit rejection cases for duplicate sample IDs and a tampered product hash. This makes failure behaviour part of the tested interface rather than an undocumented operator procedure.

## 4. FAIR and GA4GH scope

This project uses open genomic formats and machine-readable provenance in ways that support FAIR data-management principles:

- metadata and deterministic release identifiers improve findability within a controlled platform;
- standard VCF/BGEN/Parquet products improve interoperability;
- source identity, QC parameters, schemas and hashes improve reuse and reproducibility.

Access control is intentionally outside this public-data example; in a real Trusted Research Environment, accessibility would remain governed by platform policy.

The project shows awareness of GA4GH data standards but does not claim implementation of GA4GH APIs or services.
