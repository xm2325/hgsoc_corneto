# Genomic Data Engineering Pipeline

A production-style, reproducible genotype data workflow built around public 1000 Genomes Phase 3 chromosome 22 data. The pipeline validates and normalises VCF, converts genotypes to PLINK 2 PGEN, applies explicit variant QC, exports BGEN, writes analysis-ready Parquet tables, and records provenance.

## Why this exists

The project demonstrates the engineering path from an external genomic delivery to reviewable research data products:

`public VCF -> bcftools validation/normalisation -> PLINK2 PGEN -> QC -> BGEN + Parquet -> provenance`

It is intentionally a data-engineering project, not a GWAS claim. No association result or clinical conclusion is produced.

## Real public input

CI uses a public 1000 Genomes Project Phase 3 chr22 genotype subset retained in RTI International's GAWMerge test data. Its filename contains `N100`, but machine validation of the pinned VCF finds **90 sample IDs**. The pipeline therefore derives the expected sample count from the input itself and asserts that all input samples are retained downstream instead of trusting the filename. The URL is pinned to commit `ce6ad1c1af42bfa501b5e122b747f54c6644e2e9` so the bytes cannot change without an explicit source update.

Upstream call set: 1000 Genomes Phase 3, GRCh37, chromosome 22.

The first complete CI execution retained all 90 samples and produced **127,171 variants after the configured QC gates**.

## QC contract

Default filters are configuration parameters, not hidden constants:

- variant missingness: `--geno 0.02`
- minor allele frequency: `--maf 0.01`
- Hardy-Weinberg exact-test threshold: `--hwe 1e-6 midp`
- input restricted to biallelic SNPs after multiallelic decomposition

These defaults are illustrative engineering gates and are not presented as a universal scientific analysis protocol.

## Outputs

- `normalised.vcf.gz` + tabix index
- PLINK2 `qc.pgen/.pvar/.psam`
- PLINK QC reports: allele frequencies, sample and variant missingness, HWE, PCA
- Oxford BGEN 1.2 + `.sample`
- Parquet tables for variants, samples, allele frequencies, missingness, HWE and PCA
- `summary.json` with sample/variant counts and data-product inventory
- `provenance.json` with input source, SHA-256, parameters, tool versions and output hashes

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

GitHub Actions executes the full Nextflow pipeline on the pinned public genotype VCF, compares input and output sample counts, validates BGEN/Parquet/provenance products, builds the Docker image, and uploads the generated data products as a workflow artifact.

## Scope and interpretation

This repository demonstrates ingestion, data contracts, deterministic transformation, QC gates, format conversion, testability and provenance. It does not claim genotype calling, phasing, imputation, ancestry inference or clinical validity.
