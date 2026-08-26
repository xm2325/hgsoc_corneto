# Genomic Data Engineering Pipeline

A reproducible multi-feed genotype data workflow built around public 1000 Genomes Phase 3 chromosome 22 genotype data plus an independently pinned 1000 Genomes sample-metadata feed. It validates both inputs before release, normalises VCF, converts genotypes to PLINK2 PGEN, applies explicit QC, exports and re-imports BGEN under a round-trip contract, joins sample metadata under a fail-closed identity contract, writes typed query-ready Parquet, records provenance, and blocks release when source, join, interface, semantic or query contracts fail.

## Engineering path

`provider genotype manifest + VCF -> genotype delivery gate -> bcftools -> PLINK2 -> QC/PCA -> BGEN export/re-import -> BGEN contract`

`pinned sample panel -> immutable-source check -> canonical sample-ID mapping -> 100% metadata join -> metadata Parquet`

`genotype products + metadata contract -> typed Parquet -> DuckDB query contract -> provenance -> release gate -> exact resume + fresh Docker equivalence`

This is a data-engineering project, not a GWAS or clinical analysis. It does not produce an association result or clinical conclusion.

## Latest fully green real-data evidence

The latest fully green v0.6.0 code run is GitHub Actions run `32946752358` on commit `9af6e1da6f99262a4058491f357bdc57599e317d`.

- **43/43 unit and negative-path tests passed**.
- **15/15 real-data Nextflow processes passed**.
- The genotype delivery gate returned `PASS / PROCESS`; its delivery fingerprint remained `252db21a858d21a0aa0d015376d247978a355503a4d62d37421eec9061357d15`.
- The second feed is pinned to Illumina/akt commit `f1e47dd3dbbd966415b04e3b346e20fa23c29a93` and Git blob `bc447774e6bacc2f4ca3619d14bf96a1846aa4e4`.
- The metadata feed contained **2,504 rows** and matched **90/90 genotype samples (100% coverage)** using a strict canonical sample-ID contract.
- The joined metadata product had semantic SHA-256 `07190a7bf645849f2eafcfe8368eb47511387e4e2c004851ecdd9727c17d6364` and canonical ordered sample-ID SHA-256 `48b6a1ae93607949a7c1f8a96d2a06465f146a39227f64fbb2c22a7139adba1f`.
- For this pinned 90-sample fixture, joined metadata contained 73 GBR / 17 FIN samples, all EUR super-population, with 42 female / 48 male records. These are validation facts about the fixture, not analytical findings.
- **90/90 genotype samples were preserved**.
- **1,059,913 normalised variants -> 127,171 variants after configured QC**.
- Oxford BGEN 1.2 was exported with **`ref-first` allele convention and 16-bit probabilities**, then re-imported with PLINK2.
- The BGEN round trip preserved all **90 samples** and **127,171 variants**; `OBS_CT` mismatch count was 0 and maximum absolute ALT-frequency difference was **0.0** at tolerance `0.0001`.
- Seven typed, ZSTD-compressed Parquet genotype/QC tables plus one release-critical sample-metadata Parquet product were produced.
- DuckDB queried positions **16,051,249-17,051,249** and returned **1,955 variants**, matching pandas.
- **85/85 release-contract checks passed**.
- Release identity **v4** binds genotype delivery identity, QC/PCA parameters, seven genomic semantic hashes, BGEN interface semantics, and the immutable metadata source + canonical join semantics.
- Semantic release ID: `5d61489ae14365c4476795946c869578f667c96b2de9c30ff9cec7b1f424f33e`.
- An immediate identical `-resume` rerun returned **15/15 processes from cache** and left all **18 tracked release-facing outputs unchanged**.
- A fresh Docker execution passed **15/15 cross-runtime semantic-equivalence checks**, including the metadata source/join semantics, all seven genomic Parquet semantic hashes, BGEN semantics, PCA execution parameters and the same release ID.
- The uploaded evidence artifact contains **43 files / 52,844,971 bytes**, with ZIP digest `sha256:ef2e14abd4873a7d31182aa4af833b6ff9a9aa13dd4fe7c00ad0dd64eec68595`.

Numerical results are promoted here only after the corresponding branch-head GitHub Actions run is green. No production-scale throughput or clinical claim is inferred from this compact public fixture.

## Multi-feed sample metadata contract

The genotype and metadata inputs are deliberately treated as different feeds. `DOWNLOAD_SAMPLE_METADATA` downloads the panel from an immutable Git commit. `validate_sample_metadata.py` recomputes the Git blob SHA-1 and file SHA-256, requires `sample`, `pop`, `super_pop` and `gender`, rejects duplicate or empty metadata IDs, and refuses implicit ID guessing.

Because the genotype import uses PLINK2 `--double-id`, a genotype IID is accepted only when it has the exact form `<sample>_<same-sample>`. That canonical ID must appear exactly once in the metadata panel. Any malformed genotype IID, duplicate metadata ID, missing sample or coverage below 100% fails the process. The successful join preserves genotype sample order and emits `sample_metadata.parquet` plus `metadata_validation.json`.

The metadata Parquet semantic hash, canonical ordered sample hash, pinned Git blob, source SHA-256, source row count and matched-sample count are included in release identity v4. A metadata-content change therefore changes logical release identity even when genotype products are unchanged.

## Provider genotype delivery contract

The fail-closed genotype boundary in `config/delivery_manifest.json` declares a stable delivery ID, provider/source identity, expected source SHA-256, declared reference genome, expected VCF sample count and deterministic sample-roster hash. The validator recomputes source hash and sample roster from the VCF. A valid new delivery returns `PROCESS`; an exact registered duplicate returns `NOOP`; a delivery-ID collision with changed content returns `REJECT`.

## BGEN interface contract

BGEN is treated as a data interface, not merely as a file that exists. The pipeline exports Oxford BGEN 1.2 with an explicit `ref-first` allele policy and 16-bit probabilities, re-imports it with PLINK2, and compares the round-trip dataset with the source PGEN-derived dataset.

The validator checks ordered sample identity, sample hash, variant count, ordered `CHROM/POS/ID/REF/ALT` identity, frequency-table row and variant order, `OBS_CT`, and ALT-frequency drift. A failed BGEN contract blocks release.

## Query-ready data and semantic hashes

The genomic analysis layer uses stable physical types. Genomic positions and counts are `int64`; QC fractions, probabilities and PCA scores are floating point; identifiers remain strings. Files use ZSTD compression and `schema_manifest.json` records the Arrow schema for the seven genomic tables.

A canonical semantic SHA-256 is stored for every typed genomic table and separately for the joined sample-metadata table. Semantic hashes cover ordered columns, logical types, row order and canonical scalar values, but not Parquet compression bytes, row-group layout or file metadata.

## Release identity and integrity

Release identity v4 hashes the pinned genotype source, delivery fingerprint, reference genome, QC parameters, sample/variant counts, seven genomic semantic hashes, the BGEN interface contract, pinned PLINK2 PCA execution parameters (`seed=20260826`, `threads=2`, `memory=3000 MB`), and the metadata source/join contract.

Byte-level SHA-256 values for BGEN, sample, metadata and Parquet products remain in provenance and are checked for integrity, but logical identity uses semantic content where serialisation details should not matter.

## Two reproducibility contracts

- **Exact cache-resume reproducibility:** an immediate Nextflow `-resume` must return all 15 processes as `CACHED` and leave all 18 tracked release-facing file SHA-256 values unchanged.
- **Cross-runtime semantic equivalence:** a fresh Docker execution must reproduce delivery identity, counts, query result, BGEN semantics, metadata source/join semantics, seven genomic semantic hashes, PCA execution parameters and semantic release ID.

## Outputs

- genotype source VCF SHA-256 and `delivery_validation.json`;
- pinned sample metadata panel + SHA-256;
- normalised VCF + tabix index and source inventory;
- PLINK2 PGEN and QC/PCA reports;
- Oxford BGEN 1.2 + sample file + `bgen_validation.json`;
- `sample_metadata.parquet` + `metadata_validation.json`;
- seven typed ZSTD genomic Parquet tables;
- `schema_manifest.json` and `query_validation.json`;
- `summary.json`, `provenance.json`, `release_validation.json`;
- `reproducibility_validation.json` and `runtime_equivalence_validation.json`.

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

GitHub Actions executes unit and negative-path tests, the full two-feed real-data Nextflow workflow, delivery/metadata/BGEN/query/release contracts, exact cache-resume validation, Docker build, a fresh containerised real-data workflow and cross-runtime semantic-equivalence validation.

## Standards and scope

The workflow uses common genomic and analytical formats (`VCF`, `PGEN`, `BGEN`, `Parquet`) with machine-readable provenance, schemas and checksums. The design supports FAIR-style metadata and reproducibility. It does **not** claim a GA4GH API implementation, genotype calling, phasing, imputation, ancestry inference or clinical validity.
