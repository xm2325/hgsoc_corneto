#!/usr/bin/env bash
#SBATCH --job-name=reg-robustness-8
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/reg-robustness-8-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna
REG_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory
VENV_ROOT=${REPO_ROOT}/.venv-cpu
LAMBDA1=${REPO_ROOT}/data/processed/corneto/E-MTAB-14568/regulatory_robustness_8_lambda1.json
SUMMARY=${REPO_ROOT}/data/processed/corneto/E-MTAB-14568/regulatory_robustness_8_l0_vs_l1.json

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12
export PYTHONNOUSERSITE=1
export PYTHONPATH=${REPO_ROOT}/src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export GRB_LICENSE_FILE=${REPO_ROOT}/gurobi_xm.lic
cd "${REPO_ROOT}"

test ! -e "${LAMBDA1}"
test ! -e "${SUMMARY}"
srun "${VENV_ROOT}/bin/python" scripts/run_corneto_regulatory_pilot.py \
    --study E-MTAB-14568 \
    --expression "${RNA_ROOT}/aggregated/E-MTAB-14568/gene_tpm.tsv.gz" \
    --manifest data/processed/metadata/ocm_master_manifest.tsv \
    --collectri "${REG_ROOT}/collectri_tf_target.tsv.gz" \
    --pkn "${REG_ROOT}/omnipath_signed_pkn.tsv.gz" \
    --output "${LAMBDA1}" \
    --primary-only \
    --max-samples 8 \
    --min-targets 5 \
    --max-inputs 3 \
    --max-outputs 6 \
    --max-depth 3 \
    --max-edges 3000 \
    --max-flow 100 \
    --lambda-reg 1 \
    --max-seconds 300 \
    --solver gurobi

srun "${VENV_ROOT}/bin/python" scripts/compare_regulatory_robustness.py \
    --study E-MTAB-14568 \
    --receipt lambda0=data/processed/corneto/E-MTAB-14568/regulatory_pilot_bounded_gurobi_8_l0.json \
    --receipt lambda1="${LAMBDA1}" \
    --output "${SUMMARY}"
