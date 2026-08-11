#!/usr/bin/env bash
#SBATCH --job-name=reg-nmf-integration
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/reg-nmf-integration-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
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

srun "${VENV_ROOT}/bin/python" scripts/integrate_regulatory_nmf_states.py \
    --study E-MTAB-14568 \
    --regulatory data/processed/corneto/E-MTAB-14568/regulatory_full_gurobi_27_l0.json \
    --assignments data/processed/rna/E-MTAB-14568/nmf/rank_3_assignments.tsv \
    --output data/processed/corneto/E-MTAB-14568/regulatory_nmf_integration.json
