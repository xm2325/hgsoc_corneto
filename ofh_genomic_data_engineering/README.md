# Genomic Data Engineering Pipeline

A production-style, reproducible genotype data workflow built around public 1000 Genomes Phase 3 chromosome 22 data. The pipeline validates and normalises VCF, converts genotypes to PLINK 2 PGEN, applies explicit variant QC, exports BGEN, writes analysis-ready Parquet tables, records provenance, and blocks a data release when its contract fails.

## Engineering path

`public VCF -> source contract -> bcftools normalisation -> PLINK2 PGEN -> QC/PCA -> BGEN + Parquet -> provenance -> release gate`

This is a data-engineering project, not a GWAS analysis. It does not produce an association result or clinical conclusion.

## Real public input

CI uses a public 1000 Genomes Project Phase 3 chr22 genotype subset retained in RTI International's GAWMerge test data. Its filename contains `N100`, but machine inspection of the pinned VCF finds **90 sample IDs**. The workflow therefore derives its source contract from the file itself instead of trusting a filename. The upstream URL is pinned to commit `ce6ad1c1af42bfa501b5e122b747f54c6644e2e9`.

Upstream call set: 1000 Genomes Phase 3, GRCh37, chromosome 22.

The first complete CI execution retained all 90 samples and produced **127,171 variants after the configured QC gates**.

## QC contract

Default filters are configuration parameters:

- variant missingness: `--geno 0.02`
- minor allele frequency: `--maf 0.01`
- Hardy-Weinberg exact-test threshold: `--hwe 1e-6 midp`
- biallelic SNPs after multiallelic decomposition

These are example engineering gates, not a universal scientific protocol.

## Data release contract

A candidate release is accepted only when the release gate passes. It checks:

- source sample count is preserved downstream;
- released variant count is positive and cannot exceed the normalised source count;
- required Parquet tables and required columns are present;
- Parquet row counts agree with the release summary;
- sample IDs and normalised variant keys are unique;
- sample-level QC and PCA tables refer to the same sample set;
- allele frequency, missingness and HWE probability fields are bounded in `[0, 1]`;
- BGEN, sample and Parquet SHA-256 hashes agree with the provenance manifest.

A successful release receives a deterministic 64-character `release_id`. A failed contract returns a non-zero exit code, so Nextflow and GitHub Actions block publication.

See [`docs/DATA_RELEASE_CONTRACT.md`](docs/DATA_RELEASE_CONTRACT.md).

## Outputs

- `normalised.vcf.gz` + tabix index
- `source_inventory.json` with source sample/variant counts and hashes
- PLINK2 `qc.pgen/.pvar/.psam`
- PLINK QC reports: allele frequencies, sample and variant missingness, HWE, PCA
- Oxford BGEN 1.2 + `.sample`
- seven Parquet tables for variants, samples, allele frequencies, missingness, HWE and PCA
- `summary.json` with data-product row counts
- `provenance.json` with source SHA-256, QC parameters, tool versions and output hashes
- `release_validation.json` with PASS/FAIL checks and deterministic `release_id`

## Run

```bash
nextflow run main.nf -profile local
```

Containerised run:

```bash
docker build -t genomic-data-engineering .
docker run --rm -v "$PWD/results:/work/results" genomic-data-engineering
```

## Tests

```bash
python -m pip install -e '.[test]'
pytest
```

GitHub Actions executes the full Nextflow pipeline on the pinned public genotype VCF, runs positive and negative release-contract tests, validates BGEN/Parquet/provenance outputs, builds the Docker image, and uploads the generated release evidence.

## Standards and scope

The workflow uses common exchange and analysis formats (`VCF`, `PGEN`, `BGEN`, `Parquet`) and retains machine-readable provenance and checksums to support reusable data releases. The design is compatible with FAIR-style metadata and provenance principles. It does **not** claim a GA4GH API implementation, genotype calling, phasing, imputation, ancestry inference or clinical validity.
