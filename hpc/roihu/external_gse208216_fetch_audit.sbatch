#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=ext-gse208216
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-gse208216-%j.out

set -euo pipefail

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
OUT=${EXTERNAL_DATA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}
RECEIPT=${GSE208216_RECEIPT:-${OUT}/receipts/gse208216_fetch_audit.json}

cd "${ROOT}"
test -x "${PY}"
mkdir -p "${OUT}/gse208216" "${OUT}/receipts" "${ROOT}/logs"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"

srun "${PY}" scripts/prepare_organoid_depmap_validation.py fetch-gse208216 \
    --output-dir "${OUT}/gse208216" \
    --receipt "${RECEIPT}"

"${PY}" -c 'import json,sys; r=json.load(open(sys.argv[1])); assert r["scientific_success"] is True' \
    "${RECEIPT}"
