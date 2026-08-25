nextflow.enable.dsl=2

params.source_url = params.source_url ?: 'https://raw.githubusercontent.com/RTIInternational/GAWMerge/ce6ad1c1af42bfa501b5e122b747f54c6644e2e9/test_data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.N100.vcf.gz'
params.geno = params.geno ?: 0.02
params.maf  = params.maf  ?: 0.01
params.hwe  = params.hwe  ?: '1e-6'
params.outdir = params.outdir ?: 'results'

process DOWNLOAD_INPUT {
    tag '1000G chr22 N100'
    publishDir "${params.outdir}/00_source", mode: 'copy', overwrite: true
    output:
    path 'input.vcf.gz', emit: vcf
    path 'source.sha256', emit: source_sha
    script:
    """
    set -euo pipefail
    curl --fail --location --retry 3 --output input.vcf.gz '${params.source_url}'
    gzip -t input.vcf.gz
    sha256sum input.vcf.gz > source.sha256
    """
}

process NORMALISE_VCF {
    tag 'bcftools normalise'
    publishDir "${params.outdir}/01_normalised", mode: 'copy', overwrite: true
    input:
    path vcf
    output:
    path 'normalised.vcf.gz', emit: vcf
    path 'normalised.vcf.gz.tbi', emit: index
    path 'bcftools.stats.txt', emit: stats
    script:
    """
    set -euo pipefail
    bcftools view -h ${vcf} >/dev/null
    bcftools norm -m -any ${vcf} -Ou \
      | bcftools view -m2 -M2 -v snps -Oz -o normalised.vcf.gz
    bcftools index --tbi normalised.vcf.gz
    bcftools stats normalised.vcf.gz > bcftools.stats.txt
    test -s normalised.vcf.gz
    test -s normalised.vcf.gz.tbi
    """
}

process IMPORT_PLINK2 {
    tag 'VCF to PGEN'
    publishDir "${params.outdir}/02_pgen_raw", mode: 'copy', overwrite: true
    input:
    path vcf
    output:
    tuple path('raw.pgen'), path('raw.pvar'), path('raw.psam'), emit: pfile
    script:
    """
    set -euo pipefail
    plink2 --vcf ${vcf} --double-id --make-pgen --out raw
    test -s raw.pgen && test -s raw.pvar && test -s raw.psam
    """
}

process QC_FILTER {
    tag 'variant QC'
    publishDir "${params.outdir}/03_pgen_qc", mode: 'copy', overwrite: true
    input:
    tuple path(raw_pgen), path(raw_pvar), path(raw_psam)
    output:
    tuple path('qc.pgen'), path('qc.pvar'), path('qc.psam'), emit: pfile
    path 'qc.log', emit: log
    script:
    """
    set -euo pipefail
    plink2 --pfile raw \
      --geno ${params.geno} \
      --maf ${params.maf} \
      --hwe ${params.hwe} midp \
      --make-pgen --out qc
    test -s qc.pgen && test -s qc.pvar && test -s qc.psam
    """
}

process QC_REPORTS {
    tag 'QC metrics'
    publishDir "${params.outdir}/04_qc_reports", mode: 'copy', overwrite: true
    input:
    tuple path(qc_pgen), path(qc_pvar), path(qc_psam)
    output:
    path 'qc_metrics.afreq', emit: afreq
    path 'qc_metrics.vmiss', emit: vmiss
    path 'qc_metrics.smiss', emit: smiss
    path 'qc_metrics.hardy', emit: hardy
    path 'qc_metrics.eigenvec', emit: eigenvec
    path 'qc_metrics.eigenval', emit: eigenval
    path 'qc_metrics.log', emit: log
    script:
    """
    set -euo pipefail
    plink2 --pfile qc --freq --missing --hardy --pca 10 approx --out qc_metrics
    test -s qc_metrics.afreq
    test -s qc_metrics.vmiss
    test -s qc_metrics.smiss
    test -s qc_metrics.hardy
    test -s qc_metrics.eigenvec
    """
}

process EXPORT_BGEN {
    tag 'BGEN export'
    publishDir "${params.outdir}/05_bgen", mode: 'copy', overwrite: true
    input:
    tuple path(qc_pgen), path(qc_pvar), path(qc_psam)
    output:
    path 'analysis_ready.bgen', emit: bgen
    path 'analysis_ready.sample', emit: sample
    path 'analysis_ready.log', emit: log
    script:
    """
    set -euo pipefail
    plink2 --pfile qc --export bgen-1.2 --out analysis_ready
    test -s analysis_ready.bgen
    test -s analysis_ready.sample
    """
}

process EXPORT_PARQUET {
    tag 'Parquet export'
    publishDir "${params.outdir}/06_parquet", mode: 'copy', overwrite: true
    input:
    tuple path(qc_pgen), path(qc_pvar), path(qc_psam)
    path afreq
    path vmiss
    path smiss
    path hardy
    path eigenvec
    output:
    path 'parquet', emit: parquet_dir
    path 'summary.json', emit: summary
    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/metrics_to_parquet.py \
      --pvar qc.pvar --psam qc.psam \
      --afreq ${afreq} --vmiss ${vmiss} --smiss ${smiss} \
      --hardy ${hardy} --eigenvec ${eigenvec} \
      --outdir parquet --summary summary.json
    """
}

process PROVENANCE {
    tag 'provenance'
    publishDir "${params.outdir}/07_provenance", mode: 'copy', overwrite: true
    input:
    path source_sha
    path summary
    path bgen
    path sample
    path parquet_dir
    path stats
    output:
    path 'provenance.json'
    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/build_provenance.py \
      --source-url '${params.source_url}' \
      --source-sha-file ${source_sha} \
      --summary ${summary} \
      --bgen ${bgen} --sample ${sample} \
      --parquet-dir ${parquet_dir} --bcftools-stats ${stats} \
      --geno '${params.geno}' --maf '${params.maf}' --hwe '${params.hwe}' \
      --output provenance.json
    """
}

workflow {
    DOWNLOAD_INPUT()
    NORMALISE_VCF(DOWNLOAD_INPUT.out.vcf)
    IMPORT_PLINK2(NORMALISE_VCF.out.vcf)
    QC_FILTER(IMPORT_PLINK2.out.pfile)
    QC_REPORTS(QC_FILTER.out.pfile)
    EXPORT_BGEN(QC_FILTER.out.pfile)
    EXPORT_PARQUET(
        QC_FILTER.out.pfile,
        QC_REPORTS.out.afreq,
        QC_REPORTS.out.vmiss,
        QC_REPORTS.out.smiss,
        QC_REPORTS.out.hardy,
        QC_REPORTS.out.eigenvec
    )
    PROVENANCE(
        DOWNLOAD_INPUT.out.source_sha,
        EXPORT_PARQUET.out.summary,
        EXPORT_BGEN.out.bgen,
        EXPORT_BGEN.out.sample,
        EXPORT_PARQUET.out.parquet_dir,
        NORMALISE_VCF.out.stats
    )
}
