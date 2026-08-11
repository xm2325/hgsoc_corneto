#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=nmf-primary-gate

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
MATRIX_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
exec srun "$PY" scripts/check_nmf_primary_inputs.py \
  --repo "$ROOT" --matrix-root "$MATRIX_ROOT" \
  --output data/processed/rna/nmf_primary_input_gate.json
