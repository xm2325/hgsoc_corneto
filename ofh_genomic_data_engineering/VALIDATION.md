# Validation

A change is application-ready only when the same GitHub Actions job passes all branch-head gates:

1. Python unit and negative-path tests, including provider-delivery, sample-metadata join, BGEN, release and cross-runtime failure paths.
2. Installation of pinned PLINK2 and Nextflow plus runner bcftools.
3. End-to-end execution on the pinned public 1000 Genomes genotype VCF and independently pinned sample-metadata panel.
4. Genotype provider contract: source checksum, declared genome build, sample count and sample-roster hash.
5. Metadata source contract: immutable Git blob, file SHA-256, required schema and non-empty unique metadata IDs.
6. Strict canonical genotype-to-metadata identity mapping and 100% join coverage.
7. VCF normalisation, PLINK2 import and configured QC/PCA.
8. BGEN 1.2 export with explicit `ref-first` allele convention and 16-bit probabilities.
9. BGEN re-import contract covering ordered sample identity, ordered variant/allele identity, `OBS_CT` and ALT-frequency tolerance.
10. Seven typed ZSTD genomic Parquet products plus Arrow schema manifest and canonical semantic hashes.
11. Release-critical sample-metadata Parquet with canonical semantic hash.
12. DuckDB genomic region query contract cross-checked against pandas.
13. Provenance binding genotype delivery, metadata join, BGEN contract, products and pinned PLINK2 PCA execution parameters.
14. Fail-closed release contract with semantic release identity v4.
15. A second identical Nextflow run with `-resume`, where all 15 processes must be `CACHED` and all 18 tracked release-facing SHA-256 values must remain unchanged.
16. Docker image build followed by a fresh containerised run of the same two-feed real-data workflow.
17. Cross-runtime semantic-equivalence validation covering genotype delivery, metadata source/join, sample/variant counts, seven genomic semantic hashes, BGEN semantics, PCA execution parameters, query result and release ID.
18. Workflow artifact upload.

The v0.6.0 reference is GitHub Actions run `32946752358` on commit `9af6e1da6f99262a4058491f357bdc57599e317d`: **43 tests passed, 15 real-data processes passed, 85 release checks passed, 15/15 resume processes were cached, 18 tracked outputs were unchanged, and 15/15 host/Docker semantic-equivalence checks passed**.

The metadata feed contained 2,504 rows and matched 90/90 genotype samples. The joined metadata semantic SHA-256 was `07190a7bf645849f2eafcfe8368eb47511387e4e2c004851ecdd9727c17d6364`; the canonical ordered sample-ID SHA-256 was `48b6a1ae93607949a7c1f8a96d2a06465f146a39227f64fbb2c22a7139adba1f`.

Delivery negative tests verify rejection of checksum mismatch, wrong genome build, sample-roster mismatch and delivery-ID collision; exact registered duplicates return `NOOP`.

Metadata negative tests verify rejection of wrong pinned blob, missing sample, duplicate metadata ID and malformed PLINK double-ID. Release/runtime tests additionally reject metadata contract failure and metadata semantic drift.

BGEN negative tests verify rejection of sample-order drift, REF/ALT or ordered variant-identity drift, and ALT-frequency drift beyond tolerance.

Release/runtime negative tests also verify rejection of duplicate sample IDs, tampered product hashes, invalid typed schemas, declared semantic-hash mismatch, genomic semantic content drift, invalid PCA execution parameters and host/candidate PCA-parameter mismatch.

The resume check writes `results/08_release/reproducibility_validation.json`; the fresh-runtime check writes `results/08_release/runtime_equivalence_validation.json`. No README or CV numerical claim should be promoted from a branch run until the corresponding branch-head workflow is green.
