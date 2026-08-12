#!/usr/bin/env bash
# Freeze the Taylor-only regulatory signature before full external inference.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=freeze-taylor-reg
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/freeze-taylor-reg-%j.out

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
BASE=${ROOT}/data/processed/corneto/regulatory_multisample
OUT=${ROOT}/data/processed/corneto/external_validation

cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME
test ! -e "${OUT}/taylor_regulatory_signature_l0p001.tsv"
test ! -e "${OUT}/taylor_regulatory_signature_l0p001_receipt.json"
mkdir -p "${OUT}"

srun "${PY}" scripts/freeze_taylor_regulatory_signature.py \
    --pooled "${BASE}/v1/grid/pooled/l0p001.json" \
    --cohort "E-MTAB-7223=${BASE}/v1/grid/E-MTAB-7223/l0p001.json" \
    --cohort "E-MTAB-10801=${BASE}/v1/grid/E-MTAB-10801/l0p001.json" \
    --cohort "E-MTAB-11000=${BASE}/v1/grid/E-MTAB-11000/l0p001.json" \
    --cohort "E-MTAB-14568=${BASE}/v1/grid/E-MTAB-14568/l0p001.json" \
    --patient-balanced "${BASE}/patient_balanced52/regulatory_l0p001.json" \
    --richer-pooled "${BASE}/richer_v1/grid/pooled/l0p001.json" \
    --min-patient-fraction 0.10 \
    --min-cohort-fraction 0.50 \
    --output "${OUT}/taylor_regulatory_signature_l0p001.tsv" \
    --receipt "${OUT}/taylor_regulatory_signature_l0p001_receipt.json"
