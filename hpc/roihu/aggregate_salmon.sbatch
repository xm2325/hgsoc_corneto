#!/usr/bin/env bash
#SBATCH --job-name=hgsoc-aggregate
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/aggregate-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna
VENV_ROOT=${REPO_ROOT}/.venv-cpu
STUDY_ACCESSION=${STUDY_ACCESSION:?Submit with --export=ALL,STUDY_ACCESSION=E-MTAB-14568}

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12
export PYTHONNOUSERSITE=1
cd "${REPO_ROOT}"

srun "${VENV_ROOT}/bin/python" scripts/aggregate_salmon.py \
    --run-manifest data/processed/metadata/rna_runs.tsv \
    --master-manifest data/processed/metadata/ocm_master_manifest.tsv \
    --study "${STUDY_ACCESSION}" \
    --gtf "${RNA_ROOT}/reference/gencode_v32/downloads/gencode.v32.annotation.gtf.gz" \
    --quant-root "${RNA_ROOT}/quant" \
    --matrix-output-dir "${RNA_ROOT}/aggregated/${STUDY_ACCESSION}" \
    --report-output-dir "data/processed/rna/${STUDY_ACCESSION}/aggregation"
