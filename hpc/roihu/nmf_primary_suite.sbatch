#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=nmf-primary

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
MATRIX_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
case "${SLURM_ARRAY_TASK_ID:?}" in
  0) TARGET=E-MTAB-7223; OUT=data/processed/rna/E-MTAB-7223/nmf_primary ;;
  1) TARGET=E-MTAB-10801; OUT=data/processed/rna/E-MTAB-10801/nmf_primary ;;
  2) TARGET=E-MTAB-11000; OUT=data/processed/rna/E-MTAB-11000/nmf_primary ;;
  3) TARGET=E-MTAB-14568; OUT=data/processed/rna/E-MTAB-14568/nmf_primary ;;
  4) TARGET=pooled60; OUT=data/processed/rna/pooled_primary60/nmf ;;
  *) echo "unknown array index: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac
test -s data/processed/rna/nmf_primary_input_gate.json
test ! -e "$OUT"
exec srun "$PY" scripts/run_nmf_primary_suite.py \
  --repo "$ROOT" --matrix-root "$MATRIX_ROOT" --target "$TARGET" \
  --output-dir "$OUT" --top-genes 6000 --starts 100 --max-rank 10
