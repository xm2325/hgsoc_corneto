#!/usr/bin/env bash
#SBATCH --job-name=hgsoc-nmf
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/nmf-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna
VENV_ROOT=${REPO_ROOT}/.venv-cpu

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "${REPO_ROOT}"

srun "${VENV_ROOT}/bin/python" scripts/run_nmf_stability.py \
    --config config/nmf_pipeline.yaml \
    --matrix "${RNA_ROOT}/aggregated/E-MTAB-14568/gene_log1p_tpm.tsv.gz" \
    --sample-qc data/processed/rna/E-MTAB-14568/aggregation/sample_qc.tsv \
    --artifact-output-dir "${RNA_ROOT}/nmf/E-MTAB-14568/independent_k2-8_100runs" \
    --report-output-dir data/processed/rna/E-MTAB-14568/nmf
