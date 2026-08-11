#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=reg-receipt-repair

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python

ORIGINAL_RECEIPT=${ORIGINAL_RECEIPT:?}
RETRY_RECEIPT=${RETRY_RECEIPT:?}
OUTPUT_RECEIPT=${OUTPUT_RECEIPT:?}
test ! -e "$OUTPUT_RECEIPT"
exec srun "$PY" scripts/repair_regulatory_receipt.py \
  --original "$ORIGINAL_RECEIPT" \
  --retry "$RETRY_RECEIPT" \
  --output "$OUTPUT_RECEIPT"
