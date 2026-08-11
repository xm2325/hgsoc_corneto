#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=reg-single-retry

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python

STUDY=${STUDY:?}
RUN_ACCESSION=${RUN_ACCESSION:?}
LAMBDA_REG=${LAMBDA_REG:?}
OUTPUT_RECEIPT=${OUTPUT_RECEIPT:?}
EXPRESSION=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated/${STUDY}/gene_tpm.tsv.gz
SOURCE_DIR=/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory

mkdir -p "$(dirname "$OUTPUT_RECEIPT")"
test ! -e "$OUTPUT_RECEIPT"
test -s "$EXPRESSION"
exec srun "$PY" scripts/run_corneto_regulatory_pilot.py \
  --expression "$EXPRESSION" \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --collectri "$SOURCE_DIR/collectri_tf_target.tsv.gz" \
  --pkn "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz" \
  --study "$STUDY" --output "$OUTPUT_RECEIPT" --primary-only \
  --run-accession "$RUN_ACCESSION" \
  --max-inputs 3 --max-outputs 6 --max-depth 3 --max-edges 3000 --max-flow 100 \
  --lambda-reg "$LAMBDA_REG" --max-seconds 300 --solver gurobi
