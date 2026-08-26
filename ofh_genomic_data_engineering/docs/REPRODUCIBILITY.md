# Deterministic rerun and resume contract

## Goal

A successful first execution is not enough evidence for a reusable data pipeline. For a fixed source, fixed parameters and fixed code revision, an immediate rerun should not silently create a different research data release.

## CI contract

After the first real-data run passes the query and release contracts, GitHub Actions records SHA-256 values for the release-facing BGEN, sample file, seven Parquet tables, schema manifest, query validation, summary, provenance and release validation.

The workflow then executes the identical Nextflow graph with `-resume` and writes a trace file. The validation fails unless:

- the trace contains exactly the 11 expected pipeline processes;
- every process status is `CACHED`;
- every tracked output has the same SHA-256 value before and after resume;
- the deterministic release ID remains available from the release contract.

A passing check writes `results/08_release/reproducibility_validation.json`, which is uploaded with the workflow artifact.

## Interpretation

This checks deterministic immediate rerun behaviour and correct use of the Nextflow cache for the pinned CI dataset. It is not a claim that arbitrary historical pipelines remain reproducible after tool, reference-build, parameter or source changes. Those changes must create a new validated release context.
