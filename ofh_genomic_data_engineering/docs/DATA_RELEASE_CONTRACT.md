# Data Release Contract

## Purpose

The pipeline treats a genomic delivery and a research data release as two separate contracts. Processing success alone is not enough to publish a release.

## 1. Source contract

After VCF normalisation, `source_inventory.py` records sample count, an ordered sample-ID SHA-256, indexed variant count and a SHA-256 of the normalised VCF. Duplicate or zero sample IDs and zero indexed variants are rejected before release validation.

The normalised-VCF file hash is provenance/integrity evidence. It is not used as the logical release identifier because compression headers and other encoding details can differ between fresh executions.

## 2. Curated release contract

`validate_release.py` checks the BGEN, Parquet, summary and provenance products as one release. The gate validates file presence, required columns, row-count consistency, unique sample IDs, unique normalised variant keys, source-to-release sample preservation, variant-count non-inflation, bounded QC metrics, cross-table sample consistency and SHA-256 integrity.

Each of the seven Parquet tables also has a canonical semantic SHA-256. The hash is computed from ordered columns, normalised logical types and canonical scalar values in row order. It is independent of Parquet compression, row-group layout and file metadata. Release validation recomputes those hashes from the stored tables and fails if declared and observed semantic content differ.

## 3. Release identity v2

A passing release receives a 64-character release ID computed from the pinned original source SHA-256, validated delivery fingerprint, declared reference genome, QC parameters (`geno`, `maf`, `hwe`), release sample/variant counts and all seven semantic table hashes.

BGEN, sample and Parquet byte-level SHA-256 values remain in provenance and are checked for integrity, but they do not define logical identity. A regression test re-serialises unchanged table content with different Parquet compression: file SHA changes while semantic hash and release ID remain unchanged. A separate test changes logical content and requires the release ID to change.

## 4. Rejection path

The command exits non-zero when any check fails, so the Nextflow `RELEASE_GATE` fails. Negative tests cover duplicate sample IDs, tampered product hashes and declared semantic-hash mismatch.

## 5. Reproducibility scope

Two contracts are kept separate: `reproducibility_validation.json` proves exact file identity for an immediate Nextflow cache-resume execution; `runtime_equivalence_validation.json` proves logical equivalence between independent host and Docker executions by comparing delivery identity, counts, query results, semantic hashes and release ID.

## 6. FAIR and GA4GH scope

This project uses open genomic formats and machine-readable provenance in ways that support FAIR data-management principles. It does not claim implementation of GA4GH APIs or services, genotype calling, phasing, imputation, ancestry inference or clinical validity.
