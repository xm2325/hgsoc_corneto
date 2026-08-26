# Deterministic rerun and runtime-equivalence contracts

## Goal

A successful first execution is not enough evidence for a reusable multi-feed data pipeline. For fixed genotype and metadata sources, fixed parameters and fixed code, an immediate rerun should not silently create a different release. A fresh execution in a different runtime should also preserve logical release semantics even when binary serialisation details can differ.

## Exact cache-resume contract

After the first real-data run passes metadata, BGEN, query and release contracts, GitHub Actions records SHA-256 values for 18 release-facing outputs: delivery validation, sample-metadata Parquet and validation, BGEN and sample file, BGEN validation, seven genomic Parquet tables, summary, schema manifest, query validation, provenance and release validation.

The workflow reruns the identical graph with `-resume`. Validation fails unless the trace contains exactly 15 processes, every process is `CACHED`, every tracked output has the same SHA-256, and semantic release ID remains unchanged.

## Fresh host/Docker semantic-equivalence contract

GitHub Actions builds the project image and executes the complete two-feed workflow again inside a fresh container. `validate_runtime_equivalence.py` requires the same genotype delivery fingerprint, sample/variant counts, seven genomic Parquet semantic hashes, BGEN round-trip contract, metadata immutable source and join semantics, pinned PCA execution parameters, DuckDB query result, release-identity version and release ID.

Binary file hashes are intentionally not the cross-runtime equality criterion because independent serialisation can differ without logical data drift.

## v0.6.0 verified result

GitHub Actions run `32946752358` on commit `9af6e1da6f99262a4058491f357bdc57599e317d` passed both contracts:

- **15/15 processes were `CACHED`** on immediate resume;
- all **18 tracked release-facing outputs were unchanged**;
- fresh host/Docker validation passed **15/15 semantic-equivalence checks**;
- the second feed matched **90/90 samples**, and metadata semantic SHA-256 `07190a7bf645849f2eafcfe8368eb47511387e4e2c004851ecdd9727c17d6364` matched across runtimes;
- all seven genomic Parquet semantic hashes matched;
- BGEN sample and variant identity matched, with maximum absolute ALT-frequency difference `0.0`;
- release identity v4 matched with release ID `5d61489ae14365c4476795946c869578f667c96b2de9c30ff9cec7b1f424f33e`.

## Interpretation

These checks cover the pinned public CI fixture and recorded tool/parameter context. They are not a claim that arbitrary historical pipelines remain reproducible after source, tool, reference-build or parameter changes; those changes must create a new validated release context.
