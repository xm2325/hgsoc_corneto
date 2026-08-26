nextflow.enable.dsl=2

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

process DELIVERY_GATE {
    tag 'provider delivery contract'
    publishDir "${params.outdir}/00_source", mode: 'copy', overwrite: true, pattern: 'delivery_validation.json'
    input:
    path vcf
    path manifest
    output:
    path 'validated_input.vcf.gz', emit: vcf
    path 'delivery_validation.json', emit: validation
    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/validate_delivery.py \
      --manifest ${manifest} --source ${vcf} \
      --required-reference-genome '${params.required_reference_genome}' \
      --expected-source-uri '${params.source_url}' \
      --output delivery_validation.json
    python - <<'PY'
import json
payload = json.load(open('delivery_validation.json'))
assert payload['status'] == 'PASS', payload
assert payload['action'] == 'PROCESS', payload
assert payload['should_process'] is True, payload
PY
    cp ${vcf} validated_input.vcf.gz
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

process SOURCE_INVENTORY {
    tag 'source contract'
    publishDir "${params.outdir}/01_normalised", mode: 'copy', overwrite: true
    input:
    path vcf
    path index
    output:
    path 'source_inventory.json', emit: inventory
    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/source_inventory.py --vcf ${vcf} --output source_inventory.json
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
    plink2 --pfile raw --geno ${params.geno} --maf ${params.maf} --hwe ${params.hwe} midp --make-pgen --out qc
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
    plink2 --pfile qc --freq --missing --hardy --pca 10 approx \
      --seed ${params.plink_seed} --threads ${params.plink_threads} \
      --memory ${params.plink_memory_mb} require --out qc_metrics
    test -s qc_metrics.afreq && test -s qc_metrics.vmiss && test -s qc_metrics.smiss
    test -s qc_metrics.hardy && test -s qc_metrics.eigenvec
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
    plink2 --pfile qc --export bgen-1.2 ref-first bits=16 --out analysis_ready
    test -s analysis_ready.bgen && test -s analysis_ready.sample
    """
}

process BGEN_ROUNDTRIP {
    tag 'BGEN re-import contract'
    publishDir "${params.outdir}/05_bgen", mode: 'copy', overwrite: true, pattern: 'bgen_validation.json'
    input:
    tuple path(qc_pgen), path(qc_pvar), path(qc_psam)
    path source_afreq
    path bgen
    path sample
    output:
    path 'bgen_validation.json', emit: validation
    script:
    """
    set -euo pipefail
    plink2 --bgen ${bgen} ref-first --sample ${sample} --make-pgen --out bgen_roundtrip
    plink2 --pfile bgen_roundtrip --freq --out bgen_roundtrip_freq
    python ${projectDir}/scripts/validate_bgen_roundtrip.py \
      --source-pvar ${qc_pvar} --source-psam ${qc_psam} --source-afreq ${source_afreq} \
      --roundtrip-pvar bgen_roundtrip.pvar --roundtrip-psam bgen_roundtrip.psam \
      --roundtrip-afreq bgen_roundtrip_freq.afreq --frequency-tolerance 0.0001 \
      --probability-bits 16 --output bgen_validation.json
    test -s bgen_validation.json
    """
}

process EXPORT_PARQUET {
    tag 'typed Parquet export'
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
    path 'schema_manifest.json', emit: schema_manifest
    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/metrics_to_parquet.py \
      --pvar qc.pvar --psam qc.psam \
      --afreq ${afreq} --vmiss ${vmiss} --smiss ${smiss} \
      --hardy ${hardy} --eigenvec ${eigenvec} \
      --outdir parquet --summary summary.json --schema-manifest schema_manifest.json
    """
}

process QUERY_CONTRACT {
    tag 'DuckDB region query'
    publishDir "${params.outdir}/06_parquet", mode: 'copy', overwrite: true
    input:
    path parquet_dir
    path schema_manifest
    output:
    path 'query_validation.json', emit: validation
    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/query_contract.py \
      --variants ${parquet_dir}/variants.parquet \
      --schema-manifest ${schema_manifest} \
      --output query_validation.json
    """
}

process PROVENANCE {
    tag 'provenance'
    publishDir "${params.outdir}/07_provenance", mode: 'copy', overwrite: true
    input:
    path source_sha
    path delivery_validation
    path summary
    path bgen
    path sample
    path bgen_validation
    path parquet_dir
    path stats
    path schema_manifest
    path query_validation
    output:
    path 'provenance.json', emit: provenance
    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/build_provenance.py \
      --source-url '${params.source_url}' --source-sha-file ${source_sha} \
      --delivery-validation ${delivery_validation} \
      --summary ${summary} --bgen ${bgen} --sample ${sample} \
      --bgen-validation ${bgen_validation} \
      --parquet-dir ${parquet_dir} --bcftools-stats ${stats} \
      --schema-manifest ${schema_manifest} --query-validation ${query_validation} \
      --geno '${params.geno}' --maf '${params.maf}' --hwe '${params.hwe}' \
      --plink-seed '${params.plink_seed}' --plink-threads '${params.plink_threads}' \
      --plink-memory-mb '${params.plink_memory_mb}' --output provenance.json
    """
}

process RELEASE_GATE {
    tag 'release contract'
    publishDir "${params.outdir}/08_release", mode: 'copy', overwrite: true
    input:
    path source_inventory
    path delivery_validation
    path summary
    path provenance
    path bgen
    path sample
    path bgen_validation
    path parquet_dir
    path query_validation
    output:
    path 'release_validation.json', emit: validation
    script:
    """
    set -euo pipefail
    python - <<'PY'
import json
delivery = json.load(open('${delivery_validation}'))
provenance = json.load(open('${provenance}'))
query = json.load(open('${query_validation}'))
bgen_validation = json.load(open('${bgen_validation}'))
assert query['status'] == 'PASS', query
assert bgen_validation['status'] == 'PASS', bgen_validation
assert delivery['status'] == 'PASS' and delivery['action'] == 'PROCESS', delivery
assert delivery['should_process'] is True, delivery
assert provenance['source']['sha256'] == delivery['source_observed']['sha256']
assert provenance['delivery']['delivery_fingerprint'] == delivery['delivery']['delivery_fingerprint']
assert provenance['parameters']['delivery_fingerprint'] == delivery['delivery']['delivery_fingerprint']
assert provenance['delivery']['reference_genome'] == '${params.required_reference_genome}'
PY
    python ${projectDir}/scripts/validate_release.py \
      --source-inventory ${source_inventory} --summary ${summary} --provenance ${provenance} \
      --parquet-dir ${parquet_dir} --bgen ${bgen} --sample ${sample} \
      --bgen-validation ${bgen_validation} --output release_validation.json
    """
}

workflow {
    delivery_manifest = Channel.fromPath(params.delivery_manifest, checkIfExists: true)
    DOWNLOAD_INPUT()
    DELIVERY_GATE(DOWNLOAD_INPUT.out.vcf, delivery_manifest)
    NORMALISE_VCF(DELIVERY_GATE.out.vcf)
    SOURCE_INVENTORY(NORMALISE_VCF.out.vcf, NORMALISE_VCF.out.index)
    IMPORT_PLINK2(NORMALISE_VCF.out.vcf)
    QC_FILTER(IMPORT_PLINK2.out.pfile)
    QC_REPORTS(QC_FILTER.out.pfile)
    EXPORT_BGEN(QC_FILTER.out.pfile)
    BGEN_ROUNDTRIP(
        QC_FILTER.out.pfile,
        QC_REPORTS.out.afreq,
        EXPORT_BGEN.out.bgen,
        EXPORT_BGEN.out.sample
    )
    EXPORT_PARQUET(
        QC_FILTER.out.pfile,
        QC_REPORTS.out.afreq,
        QC_REPORTS.out.vmiss,
        QC_REPORTS.out.smiss,
        QC_REPORTS.out.hardy,
        QC_REPORTS.out.eigenvec
    )
    QUERY_CONTRACT(EXPORT_PARQUET.out.parquet_dir, EXPORT_PARQUET.out.schema_manifest)
    PROVENANCE(
        DOWNLOAD_INPUT.out.source_sha,
        DELIVERY_GATE.out.validation,
        EXPORT_PARQUET.out.summary,
        EXPORT_BGEN.out.bgen,
        EXPORT_BGEN.out.sample,
        BGEN_ROUNDTRIP.out.validation,
        EXPORT_PARQUET.out.parquet_dir,
        NORMALISE_VCF.out.stats,
        EXPORT_PARQUET.out.schema_manifest,
        QUERY_CONTRACT.out.validation
    )
    RELEASE_GATE(
        SOURCE_INVENTORY.out.inventory,
        DELIVERY_GATE.out.validation,
        EXPORT_PARQUET.out.summary,
        PROVENANCE.out.provenance,
        EXPORT_BGEN.out.bgen,
        EXPORT_BGEN.out.sample,
        BGEN_ROUNDTRIP.out.validation,
        EXPORT_PARQUET.out.parquet_dir,
        QUERY_CONTRACT.out.validation
    )
}
