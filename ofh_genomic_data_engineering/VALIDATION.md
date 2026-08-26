# Validation

A change is application-ready only when the same GitHub Actions job passes all branch-head gates:

1. Python unit and negative-path tests.
2. Installation of pinned PLINK2 and Nextflow plus runner bcftools.
3. End-to-end execution on the pinned public 1000 Genomes VCF.
4. Source inventory from the VCF and tabix index.
5. VCF normalisation, PLINK2 import, QC/PCA and BGEN export.
6. Seven typed ZSTD Parquet products plus Arrow schema manifest.
7. DuckDB genomic region query contract cross-checked against pandas.
8. Provenance with source/product SHA-256 values.
9. Fail-closed release contract.
10. A second identical Nextflow run with `-resume`, where all 11 processes must be `CACHED` and all tracked release-product SHA-256 values must remain unchanged.
11. Docker image build and workflow artifact upload.

Negative tests verify rejection of duplicate sample IDs, a tampered BGEN hash and an invalid string-typed genomic position schema.

The resume check writes `results/08_release/reproducibility_validation.json` with the release ID, cached-process count and post-resume product hashes. It is validation evidence rather than an input to the release ID itself.

No CV or README numerical claim should be updated from a branch run until the corresponding branch-head workflow is green.
