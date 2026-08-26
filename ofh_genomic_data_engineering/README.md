# Genomic Data Engineering Pipeline

A reproducible genotype data workflow built around public 1000 Genomes Phase 3 chromosome 22 data. It validates a provider delivery before processing, normalises VCF, converts genotypes to PLINK2 PGEN, applies explicit QC, exports and re-imports BGEN under a round-trip contract, writes typed query-ready Parquet, records provenance, and blocks a release when its data, interface or query contracts fail.

## Engineering path

`provider manifest + VCF -> delivery gate -> source contract -> bcftools -> PLINK2 -> QC/PCA -> BGEN export -> BGEN round-trip validation -> typed Parquet -> DuckDB query contract -> provenance -> release gate`

This is a data-engineering project, not a GWAS analysis. It does not produce an association result or clinical conclusion.

## Latest fully green real-data evidence

The latest fully green branch-head evidence is v0.5.0 GitHub Actions run `32944457150` on commit `bbaf6b2a6aadac5be88e4da6e65040eb00ca74d6`.

- **32/32 unit and negative-path tests passed**.
- **13/13 real-data Nextflow processes passed** on the pinned public 1000 Genomes VCF.
- The provider delivery gate returned `PASS / PROCESS` for delivery `1000g-phase3-chr22-gawmerge-n90-v1`.
- Source SHA-256, declared `GRCh37` policy, **90-sample** roster count and ordered sample-roster hash matched the manifest.
- Delivery fingerprint `252db21a858d21a0aa0d015376d247978a355503a4d62d37421eec9061357d15` was bound into provenance.
- **90/90 samples were preserved**.
- **1,059,913 normalised source variants -> 127,171 variants after configured QC**.
- BGEN 1.2 was exported with **`ref-first` allele convention and 16-bit probabilities**, then re-imported with PLINK2.
- The BGEN round-trip preserved all **90 samples** and **127,171 variants**, including ordered sample identity and ordered `CHROM/POS/ID/REF/ALT` variant identity.
- BGEN round-trip `OBS_CT` mismatches were **0** and maximum absolute ALT-frequency difference was **0.0** against a tolerance of `0.0001`.
- Seven typed, ZSTD-compressed Parquet tables plus an Arrow schema manifest were produced.
- DuckDB 1.5.5 queried positions **16,051,249-17,051,249** and returned **1,955 variants**, matching pandas.
- **71/71 release-contract checks passed**.
- Release identity **v3** binds the BGEN interface contract, semantic Parquet hashes, delivery identity, reference build, QC parameters and pinned PLINK2 PCA execution parameters.
- Semantic release ID: `556f4fb052ff3bda213d86e384e8679f3385afd8cf0ea0703d043269dd5ebd98`.
- An immediate identical `-resume` rerun returned all **13/13 processes from cache** and left all **16 tracked release-facing outputs unchanged**.
- A fresh Docker execution passed **14/14 cross-runtime semantic-equivalence checks**, including all seven Parquet semantic hashes, BGEN round-trip semantics and the same release ID.
- The workflow uploaded a **39-file, 52,832,217-byte** evidence artifact.

Numerical results are promoted here only after the corresponding branch-head GitHub Actions run is green. No latency or production-scale throughput claim is inferred from this compact public fixture.

## Provider delivery contract

The fail-closed ingestion boundary in `config/delivery_manifest.json` declares a stable delivery ID, provider/source identity, expected source SHA-256, declared reference genome, expected VCF sample count and a deterministic sample-roster hash. The validator recomputes the source hash and sample roster directly from the delivered VCF. A valid new delivery returns `PROCESS`; an exact registered duplicate returns `NOOP`; a delivery-ID collision with changed content returns `REJECT`.

## BGEN interface contract

BGEN is treated as a data interface, not merely as a file that exists. The pipeline exports Oxford BGEN 1.2 with an explicit `ref-first` allele policy and 16-bit probabilities, re-imports it with PLINK2, and compares the round-trip dataset with the source PGEN-derived dataset.

The validator checks ordered sample identity, sample hash, variant count, ordered `CHROM/POS/ID/REF/ALT` identity, frequency-table row and variant order, `OBS_CT`, and ALT-frequency drift. A passing report is written to `results/05_bgen/bgen_validation.json` and is included in provenance and release identity v3.

## Query-ready Parquet and semantic hashes

The analysis layer uses stable physical types rather than storing every PLINK field as text. Genomic positions and counts are `int64`; QC fractions, probabilities and PCA scores are floating point; identifiers remain strings. Files use ZSTD compression and `schema_manifest.json` records the Arrow schema for all seven tables.

A canonical semantic SHA-256 is stored for every typed table. The semantic hash covers ordered columns, logical types, row order and canonical scalar values, but not Parquet compression bytes, row-group layout or file metadata. Release validation recomputes every semantic hash from the stored Parquet files before accepting the candidate release.

## Data release identity and integrity

The release contract separates **file integrity** from **logical release identity**. File SHA-256 values for BGEN, sample and Parquet products detect tampering or unexpected byte changes.

Release identity v3 hashes the pinned source, delivery fingerprint, reference genome, QC parameters, sample/variant counts, seven semantic table hashes, the BGEN interface contract and pinned PLINK2 execution parameters used by approximate PCA (`seed=20260826`, `threads=2`, `memory=3000 MB`). This prevents hidden runtime state from being omitted from the logical release context.

## Two reproducibility contracts

- **Exact cache-resume reproducibility:** an immediate Nextflow `-resume` must return all 13 processes as `CACHED` and leave all 16 tracked release-facing file SHA-256 values unchanged.
- **Cross-runtime semantic equivalence:** a fresh Docker execution must reproduce the delivery fingerprint, sample/variant counts, query result, BGEN round-trip semantics, seven Parquet semantic hashes, PCA execution parameters and semantic release ID from the host run.

This keeps binary serialisation differences separate from data drift while retaining exact checksums for integrity.

## Outputs

- source VCF SHA-256 and `delivery_validation.json`;
- normalised VCF + tabix index and source inventory;
- PLINK2 PGEN and QC/PCA reports;
- Oxford BGEN 1.2 + sample file;
- `bgen_validation.json` for round-trip interface validation;
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

GitHub Actions executes unit and negative-path tests, the full real-data Nextflow workflow, delivery/BGEN/query/release contracts, exact cache-resume validation, Docker build, a fresh containerised real-data workflow and cross-runtime semantic-equivalence validation.

## Standards and scope

The workflow uses common genomic and analytical formats (`VCF`, `PGEN`, `BGEN`, `Parquet`) with machine-readable provenance, schemas and checksums. The design supports FAIR-style metadata and reproducibility. It does **not** claim a GA4GH API implementation, genotype calling, phasing, imputation, ancestry inference or clinical validity.
