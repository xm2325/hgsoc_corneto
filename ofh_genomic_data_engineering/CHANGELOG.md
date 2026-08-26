# Changelog

## 0.5.0 - 2026-08-26

- Promote BGEN from a file-presence output to a validated interface contract.
- Export Oxford BGEN 1.2 with explicit `ref-first` allele convention and 16-bit probability precision, then re-import it with PLINK2.
- Add a BGEN round-trip validator covering exact ordered sample identity, exact ordered `CHROM/POS/ID/REF/ALT` identity, aligned frequency rows, `OBS_CT` and ALT-frequency tolerance.
- Add BGEN negative tests for sample-order drift, variant/allele identity drift and frequency drift.
- Include the BGEN interface contract in provenance and release identity v3.
- Pin approximate-PCA execution conditions (`PLINK2 --seed 20260826 --threads 2 --memory 3000 require`) after a fresh Docker run exposed runtime-dependent PCA semantic drift.
- Include those PCA execution parameters in provenance, release identity and host/Docker semantic-equivalence validation.
- Add regression tests proving invalid execution parameters block release, a changed seed changes release identity, and host/candidate PCA-parameter drift fails runtime equivalence.
- Expand the real-data graph to 13 processes and exact resume evidence to 13/13 cached processes with 16 tracked release-facing outputs.
- Expand the release gate to 71 checks and host/Docker semantic equivalence to 14 checks.
- Verify v0.5.0 on real 1000 Genomes chr22 data in GitHub Actions run `32944457150`: 90 samples, 1,059,913 normalised variants, 127,171 post-QC variants, zero BGEN `OBS_CT` mismatch, zero ALT-frequency drift, seven semantic Parquet hashes matched across host and Docker, and release ID `556f4fb052ff3bda213d86e384e8679f3385afd8cf0ea0703d043269dd5ebd98`.
- Synchronise Python and Nextflow manifest versions to 0.5.0.

## 0.4.2 - 2026-08-26

- Correct the release-identity model after fresh host/container executions showed that byte-derived product hashes can change even when the logical genomic release is unchanged.
- Add canonical typed-table semantic SHA-256 values for all seven Parquet products and verify them again from the stored Parquet data at release time.
- Define release identity v2 from the pinned source, delivery fingerprint, reference genome, QC parameters, sample/variant counts and semantic table hashes; retain file SHA-256 values separately for integrity and tamper detection.
- Add negative tests for declared semantic-hash mismatch and semantic content drift, plus a regression test proving the release ID is unchanged when identical table content is serialised with different Parquet compression.
- Add a cross-runtime semantic-equivalence validator and CI evidence comparing host and Docker results independently of binary serialisation details.
- Keep the existing `-resume` gate as a separate exact-byte cache-reproducibility contract.
- Synchronise Python and Nextflow manifest versions to 0.4.2.

## 0.4.1 - 2026-08-26

- Copy the provider delivery configuration into the Docker image; v0.4.0 could build the image but did not include `config/delivery_manifest.json` required at runtime.
- Upgrade CI from Docker build-only evidence to a full containerised execution of the same pinned real-data workflow.
- Require the containerised delivery decision, delivery fingerprint, sample/variant counts and release identity to match the host workflow before CI can pass.
- Add the `python` -> `python3` runtime alias required by Nextflow processes inside the Ubuntu image.
- Synchronise Python and Nextflow manifest versions to 0.4.1.

## 0.4.0 - 2026-08-26

- Add a provider delivery manifest and a fail-closed delivery gate before genomic processing.
- Validate expected source SHA-256, declared reference genome, VCF sample count and deterministic sample-roster hash.
- Add deterministic delivery fingerprints and bind them into provenance and release parameters.
- Add idempotency semantics: an exact registered duplicate is `NOOP`; the same delivery ID with changed content is `REJECT`.
- Add six delivery-validator tests covering acceptance, checksum mismatch, wrong genome build, sample-roster mismatch, exact duplicate and delivery-ID collision.
- Expand deterministic resume validation to the 12-process graph and include delivery validation in tracked release-facing hashes.
- Synchronise Python and Nextflow manifest versions to 0.4.0.

## 0.3.1 - 2026-08-26

- Add an explicit deterministic rerun contract to GitHub Actions.
- Snapshot release-facing output SHA-256 values, rerun the same workflow with Nextflow `-resume`, and require every one of the 11 processes to be served from cache.
- Fail validation if any tracked BGEN, Parquet, schema, query, provenance or release output changes after the resumed run.
- Emit `reproducibility_validation.json` as workflow evidence and document the scope of the reproducibility claim.

## 0.3.0 - 2026-08-26

- Convert PLINK text outputs into typed Parquet instead of all-string tables.
- Store genomic positions and counts as `int64`; QC fractions, probabilities and PCA scores as floating point.
- Use ZSTD compression and emit a machine-readable Arrow schema manifest.
- Add DuckDB region-query validation cross-checked against pandas.
- Make the release stage depend on a successful query-layer contract and include schema/query evidence in provenance.

## 0.2.0 - 2026-08-26

- Add a source inventory derived from the VCF instead of provider naming.
- Add a release contract that checks schemas, row counts, sample preservation, variant-key uniqueness, QC value ranges and product hashes.
- Add release identifiers and a Nextflow release gate that exits non-zero on failed data contracts.
- Add negative tests for duplicate sample IDs and tampered product hashes.
- Document the boundary between FAIR-compatible release metadata and GA4GH API implementation.

## 0.1.1 - 2026-08-25

- Handle PLINK2 `.pvar` VCF-style metadata before the tabular header when producing Parquet data products.
- Add regression tests for metadata preambles and missing PLINK headers.
- Move Nextflow parameter defaults into `nextflow.config` to remove undefined-parameter warnings.

## 0.1.0 - 2026-08-25

- Add the initial real-data genomic engineering workflow: pinned 1000 Genomes VCF, bcftools normalisation, PLINK2 import/QC/PCA, BGEN export, Parquet outputs, provenance, Docker and GitHub Actions validation.
