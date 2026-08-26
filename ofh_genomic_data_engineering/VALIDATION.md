# Validation

A change is application-ready only when the same GitHub Actions job passes all of the following gates on the branch head:

1. Python unit and negative-path tests.
2. Installation of pinned PLINK2 and Nextflow plus the runner bcftools package.
3. End-to-end execution on the pinned public 1000 Genomes VCF.
4. Source inventory generation from the VCF itself.
5. VCF normalisation, PLINK2 import, QC/PCA, BGEN export and seven Parquet products.
6. Provenance generation with source and product SHA-256 values.
7. Release-contract PASS, including source-to-release sample preservation, schema checks, row counts, key uniqueness, bounded QC fields and product-hash verification.
8. Docker image build.
9. Workflow artifact upload.

The release gate is fail-closed: a contract failure exits non-zero and prevents a green workflow. Negative tests deliberately verify rejection of duplicate sample IDs and a tampered BGEN hash.

No CV or README numerical claim should be updated from a branch run until the corresponding branch-head workflow is green.