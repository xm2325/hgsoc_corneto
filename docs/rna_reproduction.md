# RNA-seq reproduction

## Frozen design

All four OCM studies are quantified against the same GRCh38 / GENCODE v32
reference. Salmon output is used for cross-study harmonisation, while a STAR
GENCODE-v32 count branch is retained as a method-matched sensitivity analysis
because the Tighe paper used STAR `--quantMode GeneCounts` followed by DESeq2
median-of-ratios normalization.

The independent Salmon branch is frozen in
`config/rna_pipeline.yaml`: Salmon 1.10.0, GENCODE v32 transcripts and
annotation, the GRCh38 primary assembly as full decoys, automatic library-type
detection, mapping validation, sequence-bias correction, and GC-bias
correction. The Salmon archive and all three GENCODE inputs are checksum
verified before use.

The complete public input is 444.951 GiB compressed, before indexes,
decompression, and temporary files. Do not launch it on a laptop with less than
roughly 1 TiB free working space. The workflow is intended for Slurm scratch.

## Required retained artefacts

- ENA run report and MAGE-TAB metadata snapshots;
- reference FASTA/GTF URLs and checksums;
- Salmon index metadata and version;
- one `quant.sf` and `cmd_info.json` per run;
- gene-level counts, TPM, and log1p(TPM) matrices;
- FastQC/MultiQC summaries;
- software/environment lock files.

FASTQ files may be deleted only after both ENA MD5 verification and successful
quantification have been recorded. BAM files are not a required long-term
artefact.

## Cohort construction order

1. Validate run/accession metadata.
2. Quantify every library against the same reference.
3. Aggregate transcript-to-gene matrices.
4. Apply the manifest's canonical OCM mapping.
5. Resolve duplicate passages without looking at drug response.
6. Filter to exact HGSOC tumour-only for the main cohort.
7. Infer TF activity and networks without loading Taxol labels.

## Roihu stages

The local checkout does not install Salmon or build the reference. From the
Roihu CPU login node, deploy the pinned binary and analysis dependencies once:

```bash
bash hpc/roihu/setup_rna_environment.sh
```

Build the decoy-aware reference in a scheduled job:

```bash
sbatch hpc/roihu/rna_reference.sbatch
```

After the reference receipt and job log pass validation, run the smallest
E-MTAB-14568 library as an end-to-end smoke test (array index 25,
`ERR13907062`):

```bash
sbatch --array=25 \
  --export=ALL,STUDY_ACCESSION=E-MTAB-14568 \
  hpc/roihu/salmon_quant_array.sbatch
```

Only after that run's ENA MD5 values, `quant.sf`, `meta_info.json`, detected
library type, and mapping rate are valid should the 33-run array be launched:

```bash
sbatch --array=0-32%4 \
  --export=ALL,STUDY_ACCESSION=E-MTAB-14568 \
  hpc/roihu/salmon_quant_array.sbatch
```

The `%4` concurrency cap deliberately limits simultaneous ENA transfers and
shared-filesystem pressure. Each run is idempotent: a completed receipt is
reused, while an incomplete final output is preserved for inspection rather
than silently overwritten. FASTQ transfers use checksum validation, retained
partial files, HTTP range resumption, retries for all curl transfer errors, and
a two-minute low-speed timeout so a stalled ENA connection is resumed instead
of occupying an array slot indefinitely. A full-length partial file that fails
its ENA MD5 is timestamped and quarantined, then downloaded once from byte zero;
the invalid copy remains available for audit.

Quantification requests 48 GiB per task. An initial 24 GiB request was rejected
after `ERR13907051` exceeded that limit during index loading/quantification;
verified FASTQs and the failed staging directory are retained for audit. After
the original array reaches a terminal state, collect every failed array index
and repair them in one explicitly capped array before aggregation. The recorded
E-MTAB-14568 repair used `%4`; this kept the intended four-transfer ceiling and
completed all seven missing runs. This avoids uncontrolled retries and prevents
a partly complete cohort from entering downstream stages.

Monitoring uses a fast-failure gate followed by stage-based checks. At 30--60
seconds after submission, inspect the scheduler state, the log tail, and the
first expected output marker once; this catches import, schema, path, and parser
failures before a long wait. If startup is healthy, wait for the expected stage
duration and then inspect one terminal `sacct` record together with the full log
and receipt. Do not poll the server in a tight loop.

Compact terminal-state audits for the first E-MTAB-14568 array are retained in
`data/processed/rna/etab_14568_salmon_array_summary.tsv` and
`data/processed/rna/etab_14568_salmon_array_failures.tsv`. They distinguish the
actual scheduler settings and per-task failure evidence from the recommended
submission settings above.

After all 33 runs have passed, aggregate transcripts to genes. Counts and TPM
are summed using the versioned transcript-to-gene mapping in the same frozen
GENCODE v32 GTF used to construct the Salmon reference:

GENCODE's repeated `ont` and `tag` attributes are treated as legal
multi-valued metadata. The parser retains their first value because neither is
used for transcript-to-gene assignment, while conflicting repeated identifiers
such as `gene_id` or `transcript_id` remain fatal.

```bash
sbatch --export=ALL,STUDY_ACCESSION=E-MTAB-14568 \
  hpc/roihu/aggregate_salmon.sbatch
```

The full matrices remain under
`/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated/`; the compact
sample QC and checksum receipt are written under `data/processed/rna/` for
version control. Matrix gzip streams have deterministic headers, so their
SHA-256 values are stable when the numerical content is unchanged.

## Independent E-MTAB-14568 NMF benchmark

The executable public Barnes repository does not currently expose the exact
discovery-cohort table, subtype assignments, or analysis code. The preprint's
methods state that the official analysis used DESeq2 1.26.0 VST, MAD-variable
gene lists, NMF R 0.23.0, ranks 2--10, 50 random starts for screening, and 200
starts for selected solutions. That method is recorded in
`config/nmf_pipeline.yaml` but is not claimed as reproduced here.

The frozen independent technical benchmark uses all 33 public E-MTAB-14568
runs, log1p(TPM), the 6,000 genes with greatest MAD, ranks 2--8, and 100 random
initializations per rank. Run it on Roihu after aggregation:

```bash
sbatch hpc/roihu/nmf_stability.sbatch
```

Large consensus/factor matrices stay under the Roihu RNA root. Rank metrics,
the rank-3 assignments, checksums, software versions, and the exact
configuration are retained in the repository. Deterministically ordered labels
are named `tumour_state_1`, `tumour_state_2`, and so on; they must not be renamed
Alpha/Beta/Gamma/Delta without official label alignment.

## Execution boundary

Salmon deployment, reference construction, FASTQ transfer, quantification,
aggregation, and NMF execute on Roihu. The local checkout is limited to source
editing, lightweight unit tests, review of compact receipts, and GitHub
publication. Roihu is checked only at stage boundaries (reference, one-run
smoke test, complete array, aggregation, NMF), using the scheduler record, log,
and output receipt together.
