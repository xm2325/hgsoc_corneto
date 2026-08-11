#!/bin/bash
# Response-blind per-OCM near-optimal CARNIVAL ensemble for E-MTAB-14568.
# Each array task is independently auditable and bounded to one primary OCM.
#SBATCH --job-name=reg-alt-14568
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --array=0-26%4
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/reg-alt-14568-%A_%a.out

set -euo pipefail

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
RNA_ROOT=${RNA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_rna}
SOURCE_DIR=${SOURCE_DIR:-/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
STUDY=E-MTAB-14568
LAMBDA_REG=${LAMBDA_REG:-0.1}
LABEL=${LAMBDA_LABEL:-l0p1}
OUTPUT_ROOT=${OUTPUT_ROOT:-${ROOT}/data/processed/corneto/alternative_optima/${STUDY}/${LABEL}}
OUT=${OUTPUT_ROOT}/sample_${SLURM_ARRAY_TASK_ID}.json

cd "$ROOT"
mkdir -p logs "$(dirname "$OUT")"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-${ROOT}/gurobi_xm.lic}

test -x "$PY"
test -r "$GRB_LICENSE_FILE"
test -r "${RNA_ROOT}/aggregated/${STUDY}/gene_tpm.tsv.gz"
test -r data/processed/metadata/ocm_master_manifest.tsv
test -r "${SOURCE_DIR}/collectri_tf_target.tsv.gz"
test -r "${SOURCE_DIR}/omnipath_signed_pkn.tsv.gz"
test ! -e "$OUT"

exec srun "$PY" scripts/run_regulatory_alternative_optima.py \
  --expression "${RNA_ROOT}/aggregated/${STUDY}/gene_tpm.tsv.gz" \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --collectri "${SOURCE_DIR}/collectri_tf_target.tsv.gz" \
  --pkn "${SOURCE_DIR}/omnipath_signed_pkn.tsv.gz" \
  --study "$STUDY" \
  --output "$OUT" \
  --primary-only \
  --sample-index "$SLURM_ARRAY_TASK_ID" \
  --lambda-reg "$LAMBDA_REG" \
  --percentage 0.03 \
  --scale 0.03 \
  --rel-opt-tol 0.05 \
  --max-alternatives 30 \
  --max-seconds 300 \
  --solver gurobi
