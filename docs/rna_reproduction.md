# RNA-seq reproduction

## Frozen design

All four OCM studies are quantified against the same GRCh38 / GENCODE v32
reference. Salmon output is used for cross-study harmonisation, while a STAR
GENCODE-v32 count branch is retained as a method-matched sensitivity analysis
because the Tighe paper used STAR `--quantMode GeneCounts` followed by DESeq2
median-of-ratios normalization.

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
