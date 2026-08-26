# Deterministic rerun and runtime-equivalence contracts

## Goal

A successful first execution is not enough evidence for a reusable data pipeline. For a fixed source, fixed parameters and fixed code revision, an immediate rerun should not silently create a different research data release. A fresh execution in a different runtime should also preserve the logical release even when binary serialisation details are allowed to differ.

## Exact cache-resume contract

After the first real-data run passes the BGEN, query and release contracts, GitHub Actions records SHA-256 values for 16 release-facing outputs: delivery validation, BGEN and sample file, BGEN validation, seven Parquet tables, schema manifest, query validation, summary, provenance and release validation.

The workflow then executes the identical Nextflow graph with `-resume` and writes a trace file. Validation fails unless:

- the trace contains exactly the 13 expected pipeline processes;
- every process status is `CACHED`;
- every tracked output has the same SHA-256 value before and after resume;
- the semantic release ID remains unchanged.

A passing check writes `results/08_release/reproducibility_validation.json`.

## Fresh host/Docker semantic-equivalence contract

GitHub Actions also builds the project image and executes the complete real-data workflow again inside a fresh container. `validate_runtime_equivalence.py` compares the independent host and container results.

The check requires the same delivery fingerprint, sample and variant counts, seven Parquet semantic hashes, BGEN round-trip contract, pinned PCA execution parameters, DuckDB query result, release-identity version and semantic release ID.

Binary file hashes are intentionally not the cross-runtime equality criterion because independent Parquet/BGEN serialisation can differ without changing logical content.

## v0.5.0 verified result

GitHub Actions run `32944457150` passed both contracts:

- `13/13` processes were `CACHED` on immediate resume;
- all `16` tracked release-facing outputs were unchanged;
- fresh host/Docker validation passed `14/14` semantic-equivalence checks;
- all seven Parquet semantic hashes matched;
- BGEN sample and variant identity matched;
- BGEN maximum absolute ALT-frequency difference was `0.0`;
- release identity v3 matched with release ID `556f4fb052ff3bda213d86e384e8679f3385afd8cf0ea0703d043269dd5ebd98`.

## Interpretation

These checks cover the pinned public CI dataset and recorded tool/parameter context. They are not a claim that arbitrary historical pipelines remain reproducible after tool, reference-build, parameter or source changes. Those changes must create a new validated release context.
