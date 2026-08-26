# Genomic Data Engineering Pipeline

A reproducible genotype data workflow built around public 1000 Genomes Phase 3 chromosome 22 data. It validates a provider delivery before processing, normalises VCF, converts genotypes to PLINK2 PGEN, applies explicit QC, exports BGEN, writes typed query-ready Parquet, records provenance, and blocks a release when its data or query contracts fail.

## Engineering path

`provider manifest + VCF -> delivery gate -> source contract -> bcftools -> PLINK2 -> QC/PCA -> BGEN + typed Parquet -> DuckDB query contract -> provenance -> release gate`

This is a data-engineering project, not a GWAS analysis. It does not produce an association result or clinical conclusion.

## Latest fully green real-data evidence

The latest fully green published evidence remains the v0.4.0 branch-head run `32939354889` while v0.4.2 is being validated. That run established the real-data baseline:

- **16/16 unit and negative/idempotency tests passed**;
- **12/12 real-data Nextflow processes passed**, including the provider delivery gate;
- the real delivery returned `PASS / PROCESS` for delivery `1000g-phase3-chr22-gawmerge-n90-v1`;
- source SHA-256, declared `GRCh37` policy, **90-sample** roster count and sample-roster hash matched the manifest;
- delivery fingerprint `252db21a858d21a0aa0d015376d247978a355503a4d62d37421eec9061357d15` was bound into provenance;
- **90/90 samples were preserved**;
- **1,059,913 normalised source variants -> 127,171 variants after configured QC**;
- seven typed, ZSTD-compressed Parquet tables plus an Arrow schema manifest were produced;
- DuckDB 1.5.5 queried positions **16,051,249-17,051,249** and returned **1,955 variants**, matching pandas;
- **48 release-contract checks passed**;
- an immediate identical `-resume` rerun returned all **12/12 processes from cache** and left all **15 tracked release-facing SHA-256 values unchanged**;
- BGEN 1.2, provenance, Docker image build and workflow artifact upload passed.

The byte-derived release ID from that older run is intentionally not promoted as a stable cross-execution identifier. Fresh executions later showed that binary encoding can change without changing the logical release. v0.4.2 corrects that identity model.

Numerical results are only promoted here after the corresponding branch-head GitHub Actions run is green. No latency or production-scale throughput claim is inferred from this compact public fixture.

## Provider delivery contract

The fail-closed ingestion boundary in `config/delivery_manifest.json` declares a stable delivery ID, provider/source identity, expected source SHA-256, declared reference genome, expected VCF sample count and a deterministic sample-roster hash. The validator recomputes the source hash and sample roster directly from the delivered VCF. A valid new delivery returns `PROCESS`; an exact registered duplicate returns `NOOP`; a delivery-ID collision with changed content returns `REJECT`.

## Query-ready Parquet and semantic hashes

The analysis layer uses stable physical types rather than storing every PLINK field as text. Genomic positions and counts are `int64`; QC fractions, probabilities and PCA scores are floating point; identifiers remain strings. Files use ZSTD compression and `schema_manifest.json` records the Arrow schema for all seven tables.

v0.4.2 adds a canonical semantic SHA-256 for every typed table. The semantic hash covers ordered columns, logical types, row order and canonical scalar values, but not Parquet compression bytes, row-group layout or file metadata. Release validation recomputes every semantic hash from the stored Parquet files before accepting the candidate release.

## Data release identity and integrity

The release contract deliberately separates **file integrity** from **logical release identity**. File SHA-256 values for BGEN, sample and Parquet products detect tampering or unexpected byte changes. Release identity v2 instead hashes the pinned source, delivery fingerprint, reference genome, QC parameters, sample/variant counts and seven semantic table hashes.

Changing logical table content changes the release ID, while re-serialising identical logical table content with different Parquet compression does not. File hashes still change in the latter case and remain visible in provenance.

## Two reproducibility contracts

- **Exact cache-resume reproducibility:** an immediate Nextflow `-resume` must return all 12 processes as `CACHED` and leave tracked release-facing file SHA-256 values unchanged.
- **Cross-runtime semantic equivalence:** a fresh Docker execution must reproduce the delivery fingerprint, sample/variant counts, query result, seven semantic hashes and semantic release ID from the host run.

This prevents binary serialisation details from being mistaken for data drift while retaining exact checksums for integrity.

## Outputs

- source VCF SHA-256 and `delivery_validation.json`;
- normalised VCF + tabix index and source inventory;
- PLINK2 PGEN and QC/PCA reports;
- Oxford BGEN 1.2 + sample file;
- seven typed ZSTD Parquet tables;
- `schema_manifest.json` and `query_validation.json`;
- `summary.json`, `provenance.json`, `release_validation.json`;
- `reproducibility_validation.json` for exact cache-resume evidence;
- `runtime_equivalence_validation.json` for host/Docker semantic equivalence.

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

GitHub Actions executes unit and negative-path tests, the full real-data Nextflow workflow, delivery/query/release contracts, exact cache-resume validation, Docker build, a fresh containerised real-data workflow and cross-runtime semantic-equivalence validation.

## Standards and scope

The workflow uses common genomic and analytical formats (`VCF`, `PGEN`, `BGEN`, `Parquet`) with machine-readable provenance, schemas and checksums. The design supports FAIR-style metadata and reproducibility. It does **not** claim a GA4GH API implementation, genotype calling, phasing, imputation, ancestry inference or clinical validity.
