#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --job-name=ext-depmap-gate
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-depmap-gate-%j.out

set -euo pipefail

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
OUT=${EXTERNAL_DATA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}
RECEIPT=${DEPMAP_PREFLIGHT_RECEIPT:-${OUT}/receipts/depmap_preflight.json}

cd "${ROOT}"
test -x "${PY}"
mkdir -p "${OUT}/receipts" "${ROOT}/logs"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"

args=(
    depmap-preflight
    --receipt "${RECEIPT}"
)
if [[ -n "${DEPMAP_RELEASE:-}" ]]; then
    args+=(--release "${DEPMAP_RELEASE}")
fi
if [[ -n "${DEPMAP_MODEL_CSV:-}" ]]; then
    args+=(--models "${DEPMAP_MODEL_CSV}")
fi
if [[ -n "${DEPMAP_GENE_EFFECT_CSV:-}" ]]; then
    args+=(--gene-effect "${DEPMAP_GENE_EFFECT_CSV}")
fi
if [[ -n "${DEPMAP_RELEASE_README:-}" ]]; then
    args+=(--release-readme "${DEPMAP_RELEASE_README}")
fi

srun "${PY}" scripts/prepare_organoid_depmap_validation.py "${args[@]}"
