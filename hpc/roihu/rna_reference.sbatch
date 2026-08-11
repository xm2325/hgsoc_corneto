#!/usr/bin/env bash
#SBATCH --job-name=hgsoc-rna-reference
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/rna-reference-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
VENV_ROOT=${REPO_ROOT}/.venv-cpu
SALMON=/scratch/project_2012997/xiaomei/hgsoc_corneto_tools/salmon-1.10.0/bin/salmon

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

srun "${VENV_ROOT}/bin/python" scripts/prepare_salmon_reference.py \
    --config config/rna_pipeline.yaml \
    --salmon "${SALMON}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --receipt data/processed/rna/gencode_v32_salmon_index_receipt.json
