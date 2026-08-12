#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=ext-gse189-cpm
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-gse189-cpm-%j.out

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
BASE=${EXTERNAL_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}/GSE189955/processed
cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME

srun "${PY}" scripts/normalize_external_pseudobulk_cpm.py \
    --input "${BASE}/patient_celltype_site_raw_counts.tsv.gz" \
    --output "${BASE}/patient_celltype_site_cpm.tsv.gz" \
    --receipt "${BASE}/patient_celltype_site_cpm_receipt.json" \
    --dataset GSE189955 --expected-samples 120 --expected-genes 24153
