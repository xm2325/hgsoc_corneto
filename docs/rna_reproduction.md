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
than silently overwritten.

Monitoring remains stage-based: record the job ID, wait for the expected stage
duration, then inspect one terminal `sacct` record together with the job log and
receipt. Do not poll the server in a tight loop.
