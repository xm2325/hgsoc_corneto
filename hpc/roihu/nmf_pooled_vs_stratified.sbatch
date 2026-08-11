#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=nmf-pool-v-strat

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
exec srun "$PY" scripts/compare_nmf_pooled_stratified.py \
  --repo "$ROOT" --output-dir data/processed/rna/nmf_pooled_vs_stratified
