#!/bin/bash
# One-sample, response-blind API smoke for CORNETO near-optimal sampling.
# This is intentionally bounded; it is not the production ensemble run.
#SBATCH --job-name=reg-alt-smoke
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/reg-alt-smoke-%j.out

set -euo pipefail

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
RNA_ROOT=${RNA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_rna}
SOURCE_DIR=${SOURCE_DIR:-/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
STUDY=${STUDY:-E-MTAB-14568}
OUT=${OUTPUT:-${ROOT}/data/processed/corneto/smoke/${STUDY}_alternative_optima.json}

cd "$ROOT"
mkdir -p logs "$(dirname "$OUT")"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-${ROOT}/gurobi_xm.lic}

test -x "$PY"
test -r "$GRB_LICENSE_FILE"
test -r "${RNA_ROOT}/aggregated/${STUDY}/gene_tpm.tsv.gz"
test -r data/processed/metadata/ocm_master_manifest.tsv
test -r "${SOURCE_DIR}/collectri_tf_target.tsv.gz"
test -r "${SOURCE_DIR}/omnipath_signed_pkn.tsv.gz"

exec srun "$PY" scripts/run_regulatory_alternative_optima.py \
  --expression "${RNA_ROOT}/aggregated/${STUDY}/gene_tpm.tsv.gz" \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --collectri "${SOURCE_DIR}/collectri_tf_target.tsv.gz" \
  --pkn "${SOURCE_DIR}/omnipath_signed_pkn.tsv.gz" \
  --study "$STUDY" \
  --output "$OUT" \
  --primary-only \
  --max-samples 1 \
  --lambda-reg 0.1 \
  --percentage 0.03 \
  --scale 0.03 \
  --rel-opt-tol 0.05 \
  --max-alternatives 4 \
  --max-seconds 120 \
  --solver gurobi
