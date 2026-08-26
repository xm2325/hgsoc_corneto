# Data Release Contract

## Purpose

The pipeline treats genotype delivery, independent sample metadata, format/interface transformation and the final research-data release as separate contracts. Process completion alone is not enough to publish a release.

## 1. Genotype source contract

The provider genotype manifest declares stable delivery identity, source SHA-256, declared genome build and expected sample roster. The delivery validator recomputes source identity before processing. After normalisation, `source_inventory.py` records sample count, ordered sample-ID SHA-256, indexed variant count and normalised-VCF SHA-256.

## 2. Sample metadata source and join contract

The second feed is pinned to an immutable Git commit and expected Git blob SHA-1. `validate_sample_metadata.py` verifies the blob and file SHA-256, required fields (`sample`, `pop`, `super_pop`, `gender`), non-empty unique metadata IDs and the genotype sample IDs supplied by PLINK2.

Because genotype import uses `--double-id`, the only accepted mapping is an IID of exact form `<sample>_<same-sample>` to canonical `<sample>`. The validator rejects malformed IDs, duplicate metadata IDs and any missing genotype sample. Release requires 100% coverage and preserved genotype order.

A passing join emits `sample_metadata.parquet` and `metadata_validation.json`. The latter records the immutable source identity, 2,504-row source count for the pinned fixture, 90/90 coverage, canonical ordered sample-ID SHA-256 and a canonical semantic SHA-256 of the joined output.

## 3. BGEN round-trip contract

The pipeline exports Oxford BGEN 1.2 from the QC PGEN dataset with explicit `ref-first` allele convention and 16-bit probability precision, re-imports it with PLINK2, and checks exact ordered sample and `CHROM/POS/ID/REF/ALT` identity, aligned frequency rows, `OBS_CT`, and ALT-frequency tolerance.

A failed BGEN contract blocks release.

## 4. Curated release contract

`validate_release.py` checks BGEN, metadata, Parquet, summary and provenance products as one release. The gate validates file presence, required columns, row counts, unique sample IDs and variant keys, source-to-release sample preservation, variant-count non-inflation, bounded QC metrics, cross-table sample consistency, metadata 100% coverage and ordered sample identity, BGEN contract status, semantic hashes, provenance bindings and byte-level product integrity.

## 5. Semantic release identity v4

A passing release receives a 64-character release ID computed from:

- pinned original genotype source SHA-256;
- validated genotype delivery fingerprint;
- declared reference genome;
- QC parameters (`geno`, `maf`, `hwe`);
- pinned PLINK2 approximate-PCA execution conditions: seed, thread count and memory budget;
- release sample/variant counts;
- all seven genomic semantic table hashes;
- BGEN contract identity: format, allele convention, probability bits, frequency tolerance, ordered sample hash and ordered variant/allele hash;
- metadata contract identity: contract version, immutable source Git blob SHA-1, metadata source SHA-256, source row count, matched-sample count, canonical ordered sample hash and joined-output semantic hash.

BGEN, sample, metadata and Parquet byte-level SHA-256 values remain in provenance and are checked for integrity, but they do not define logical identity where serialisation can vary without changing data meaning.

For v0.6.0 real-data run `32946752358`, all **85 release checks passed** and release identity v4 was `5d61489ae14365c4476795946c869578f667c96b2de9c30ff9cec7b1f424f33e`.

## 6. Rejection path

The command exits non-zero when any check fails, so the Nextflow `RELEASE_GATE` fails. Negative tests cover genotype-delivery drift, malformed or incomplete metadata joins, metadata semantic drift, BGEN drift, duplicate sample IDs, tampered products, declared genomic semantic-hash mismatch, genomic semantic content drift and invalid PCA execution parameters.

## 7. Reproducibility scope

Two contracts stay separate. `reproducibility_validation.json` proves exact file identity for an immediate Nextflow cache-resume execution; v0.6.0 requires 15 cached processes and 18 unchanged tracked outputs. `runtime_equivalence_validation.json` proves logical equivalence between independent host and Docker executions by comparing genotype delivery identity, metadata source/join semantics, counts, query results, BGEN semantics, PCA execution parameters, genomic semantic hashes and release ID. v0.6.0 passed all 15 cross-runtime checks.

## 8. FAIR and GA4GH scope

This project uses open genomic formats and machine-readable provenance in ways that support FAIR data-management principles. It does not claim implementation of GA4GH APIs or services, genotype calling, phasing, imputation, ancestry inference or clinical validity.
