# Data contract

## Input

A gzip-compressed, genotype-bearing VCF readable by bcftools and PLINK2. Sample identifiers must be present in the VCF header. The pipeline records a SHA-256 digest before transformation.

## Normalised genotype layer

The current workflow decomposes multiallelic records and retains biallelic SNPs. It deliberately does not perform reference-based left alignment because the pinned CI fixture does not ship a reference FASTA. A production deployment should provide and checksum the matching reference assembly before enabling `bcftools norm -f`.

## QC layer

QC thresholds are command-line parameters. The PGEN output must be non-empty after filtering. PLINK2 reports are retained as first-class data products rather than parsed and discarded.

## Research-ready exports

BGEN 1.2 is emitted for statistical-genetics interoperability. Parquet exports preserve the PLINK report tables as typed columnar datasets for data-lake or Spark/Databricks use.

## Provenance

The manifest includes source URL and digest, QC parameters, tool version strings, summary counts, output byte sizes and SHA-256 hashes. A downstream consumer can therefore distinguish data changes, configuration changes and software changes.
