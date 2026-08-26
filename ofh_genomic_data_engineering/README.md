# Genomic Data Engineering Pipeline

A reproducible genotype data workflow built around public 1000 Genomes Phase 3 chromosome 22 data. It validates and normalises VCF, converts genotypes to PLINK2 PGEN, applies explicit QC, exports BGEN, writes typed query-ready Parquet, records provenance, and blocks a release when its data or query contracts fail.

## Engineering path

`public VCF -> source contract -> bcftools -> PLINK2 -> QC/PCA -> BGEN + typed Parquet -> DuckDB query contract -> provenance -> release gate`

This is a data-engineering project, not a GWAS analysis. It does not produce an association result or clinical conclusion.

## Verified real-data evidence (v0.3.0)

The public source is a pinned 1000 Genomes Project Phase 3 chr22 genotype subset retained in RTI International's GAWMerge test data. Its filename contains `N100`, but direct inspection finds **90 sample IDs**; the workflow therefore derives the source contract from the VCF itself.

The v0.3.0 branch-head GitHub Actions run validated:

- **10/10 unit and negative-path tests passed**;
- **11/11 real-data Nextflow processes passed**;
- **90/90 samples preserved**;
- **1,059,913 normalised source variants -> 127,171 variants after configured QC**;
- seven typed, ZSTD-compressed Parquet tables plus an Arrow schema manifest;
- DuckDB 1.5.5 genomic region query over positions **16,051,249–17,051,249** returned **1,955 variants**, with the count independently matched by pandas;
- **48/48 release-contract checks passed**;
- BGEN 1.2, provenance, Docker image build and workflow artifact upload all passed.

Numerical results are only added here after the corresponding branch-head GitHub Actions run is green. No latency or large-scale throughput claim is inferred from this compact public fixture.

## QC contract

Default configurable filters are variant missingness `--geno 0.02`, minor allele frequency `--maf 0.01`, Hardy-Weinberg exact-test threshold `--hwe 1e-6 midp`, and biallelic SNPs after multiallelic decomposition. These are example engineering gates, not a universal scientific protocol.

## Query-ready Parquet

The analysis layer uses stable physical types rather than storing every PLINK field as text. Genomic positions and counts are `int64`; QC fractions, probabilities and PCA scores are floating point; identifiers remain strings. Files use ZSTD compression and `schema_manifest.json` records the Arrow schema for all seven tables.

A DuckDB process executes a genomic range predicate against `variants.parquet` and cross-checks the count against pandas. Schema validation is fail-closed: an invalid genomic position type returns a structured contract failure before SQL execution. Any type or query-contract failure stops the pipeline before a release can pass. See [`docs/QUERY_READY_DATA.md`](docs/QUERY_READY_DATA.md).

## Data release contract

A candidate release is accepted only when the release gate passes checks for source sample preservation, variant-count non-inflation, required tables/columns, row-count consistency, key uniqueness, bounded QC metrics, cross-table sample consistency and SHA-256 integrity. It receives a deterministic 64-character `release_id` only after those checks pass.

Negative-path tests explicitly verify rejection of duplicate sample IDs, a tampered BGEN hash and an invalid string-typed genomic position schema. See [`docs/DATA_RELEASE_CONTRACT.md`](docs/DATA_RELEASE_CONTRACT.md).

## Outputs

- normalised VCF + tabix index and source inventory;
- PLINK2 PGEN and QC/PCA reports;
- Oxford BGEN 1.2 + sample file;
- seven typed ZSTD Parquet tables;
- `schema_manifest.json` and `query_validation.json`;
- `summary.json`, `provenance.json` and `release_validation.json`.

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

GitHub Actions executes unit and negative-path tests, the full real-data Nextflow workflow, query/release contracts, Docker build and artifact upload.

## Standards and scope

The workflow uses common genomic and analytical formats (`VCF`, `PGEN`, `BGEN`, `Parquet`) with machine-readable provenance, schemas and checksums. The design supports FAIR-style metadata and reproducibility. It does **not** claim a GA4GH API implementation, genotype calling, phasing, imputation, ancestry inference or clinical validity.
