#!/usr/bin/env bash
#SBATCH --job-name=hgsoc-salmon
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/salmon-%A_%a.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna
VENV_ROOT=${REPO_ROOT}/.venv-cpu
SALMON=/scratch/project_2012997/xiaomei/hgsoc_corneto_tools/salmon-1.10.0/bin/salmon
STUDY_ACCESSION=${STUDY_ACCESSION:?Submit with --export=ALL,STUDY_ACCESSION=E-MTAB-14568}

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "${REPO_ROOT}"

srun "${VENV_ROOT}/bin/python" scripts/run_salmon_quant.py \
    --manifest data/processed/metadata/rna_runs.tsv \
    --study "${STUDY_ACCESSION}" \
    --array-index "${SLURM_ARRAY_TASK_ID}" \
    --salmon "${SALMON}" \
    --index "${RNA_ROOT}/reference/gencode_v32/salmon_index" \
    --fastq-root "${RNA_ROOT}/fastq" \
    --quant-root "${RNA_ROOT}/quant" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --library-type A
