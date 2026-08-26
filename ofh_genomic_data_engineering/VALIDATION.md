# Validation

A change is application-ready only when the same GitHub Actions job passes all branch-head gates:

1. Python unit and negative-path tests, including provider-delivery rejection/idempotency paths.
2. Installation of pinned PLINK2 and Nextflow plus runner bcftools.
3. End-to-end execution on the pinned public 1000 Genomes VCF.
4. Provider delivery contract: source checksum, declared genome build, sample count and sample-roster hash.
5. Source inventory from the normalised VCF and tabix index.
6. VCF normalisation, PLINK2 import, QC/PCA and BGEN export.
7. Seven typed ZSTD Parquet products plus Arrow schema manifest.
8. DuckDB genomic region query contract cross-checked against pandas.
9. Provenance with source/product SHA-256 values and validated delivery fingerprint.
10. Fail-closed release contract with delivery/provenance consistency assertions.
11. A second identical Nextflow run with `-resume`, where all 12 processes must be `CACHED` and all tracked release-facing SHA-256 values must remain unchanged.
12. Docker image build and workflow artifact upload.

Delivery negative tests verify rejection of checksum mismatch, wrong declared genome build, sample-roster mismatch and delivery-ID collision. They also verify that an exact registered duplicate returns `NOOP` with `should_process=false`.

Existing data-product negative tests verify rejection of duplicate sample IDs, a tampered BGEN hash and an invalid string-typed genomic position schema.

The resume check writes `results/08_release/reproducibility_validation.json` with the release ID, delivery fingerprint, cached-process count and post-resume product hashes. It is validation evidence rather than an input to the release ID itself.

No CV or README numerical claim should be updated from a branch run until the corresponding branch-head workflow is green.
