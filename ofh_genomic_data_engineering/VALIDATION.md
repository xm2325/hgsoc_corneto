# Validation

A change is application-ready only when the same GitHub Actions job passes all branch-head gates:

1. Python unit and negative-path tests, including provider-delivery rejection/idempotency and BGEN round-trip failure paths.
2. Installation of pinned PLINK2 and Nextflow plus runner bcftools.
3. End-to-end execution on the pinned public 1000 Genomes VCF.
4. Provider delivery contract: source checksum, declared genome build, sample count and sample-roster hash.
5. Source inventory from the normalised VCF and tabix index.
6. VCF normalisation, PLINK2 import and configured QC/PCA.
7. BGEN 1.2 export with explicit `ref-first` allele convention and 16-bit probabilities.
8. BGEN re-import and round-trip contract covering ordered sample identity, ordered variant/allele identity, `OBS_CT` and ALT-frequency tolerance.
9. Seven typed ZSTD Parquet products plus Arrow schema manifest and canonical semantic hashes.
10. DuckDB genomic region query contract cross-checked against pandas.
11. Provenance with source/product SHA-256 values, validated delivery fingerprint, BGEN contract and pinned PLINK2 PCA execution parameters.
12. Fail-closed release contract with delivery/provenance/BGEN consistency assertions and semantic release identity v3.
13. A second identical Nextflow run with `-resume`, where all 13 processes must be `CACHED` and all 16 tracked release-facing SHA-256 values must remain unchanged.
14. Docker image build followed by a fresh containerised run of the same real-data workflow.
15. Cross-runtime semantic-equivalence validation covering delivery identity, sample/variant counts, seven semantic hashes, BGEN round-trip semantics, PCA execution parameters, query result and release ID.
16. Workflow artifact upload.

The v0.5.0 branch-head reference is GitHub Actions run `32944457150`: 32 tests passed, 13 real-data processes passed, 71 release checks passed, 13/13 resume processes were cached, and 14/14 host/Docker semantic-equivalence checks passed.

Delivery negative tests verify rejection of checksum mismatch, wrong declared genome build, sample-roster mismatch and delivery-ID collision. They also verify that an exact registered duplicate returns `NOOP` with `should_process=false`.

BGEN negative tests verify rejection of sample-order drift, REF/ALT or ordered variant-identity drift, and ALT-frequency drift beyond the declared tolerance.

Release/runtime negative tests verify rejection of duplicate sample IDs, tampered product hashes, invalid typed schemas, declared semantic-hash mismatch, semantic content drift, invalid PCA execution parameters and host/candidate PCA-parameter mismatch.

The resume check writes `results/08_release/reproducibility_validation.json` with the release ID, delivery fingerprint, BGEN variant identity, cached-process count and post-resume product hashes. The fresh-runtime check writes `results/08_release/runtime_equivalence_validation.json`.

No CV or README numerical claim should be updated from a branch run until the corresponding branch-head workflow is green.
