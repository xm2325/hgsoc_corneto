#!/usr/bin/env bash
#SBATCH --job-name=ext-gse189955
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-gse189955-%j.out

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
EXTERNAL_ROOT=${EXTERNAL_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}
PYTHON_BIN=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
FINAL_DIR=${EXTERNAL_ROOT}/GSE189955
STAGE_DIR=${EXTERNAL_ROOT}/.staging/GSE189955.${SLURM_JOB_ID:?}

cd "${REPO_ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${REPO_ROOT}/src"

# This job only downloads and aggregates public RNA counts.  It must not open
# a Gurobi session or depend on a solver license.
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME

test ! -e "${FINAL_DIR}"
test ! -e "${STAGE_DIR}"
mkdir -p "${EXTERNAL_ROOT}/.staging"
mkdir "${STAGE_DIR}"

srun "${PYTHON_BIN}" scripts/prepare_gse189955_external_baseline.py \
    --fetch \
    --raw-dir "${STAGE_DIR}/raw" \
    --output-dir "${STAGE_DIR}/processed"

test -s "${STAGE_DIR}/processed/receipt.json"
test -s "${STAGE_DIR}/processed/corneto_manifest.tsv"
test -s "${STAGE_DIR}/processed/patient_celltype_site_raw_counts.tsv.gz"

"${PYTHON_BIN}" -c \
    'import json, pathlib, sys; p=pathlib.Path(sys.argv[1]); r=json.loads(p.read_text()); assert r["status"] == "valid"; assert r["corneto_manifest"]["rows"] == r["aggregation"]["groups"]; assert r["corneto_manifest"]["definitive_malignant_rows"] == 0' \
    "${STAGE_DIR}/processed/receipt.json"

{
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID}"
    printf 'repo_commit\t%s\n' "$(git rev-parse HEAD)"
    printf 'published_path\t%s\n' "${FINAL_DIR}"
} > "${STAGE_DIR}/execution.tsv"

# STAGE_DIR and FINAL_DIR share EXTERNAL_ROOT, so this rename is atomic.  A
# failed job deliberately leaves its uniquely named stage directory for audit.
mv "${STAGE_DIR}" "${FINAL_DIR}"
