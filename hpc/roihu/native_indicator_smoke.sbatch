#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --job-name=native-link-smoke
set -euo pipefail
ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
cd "$ROOT"
export PYTHONPATH="$ROOT/scripts:$ROOT/src" PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
export GRB_LICENSE_FILE="${ROOT}/gurobi_xm.lic"
test -r "$GRB_LICENSE_FILE"
exec srun /scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python \
 scripts/validate_binary_support.py --root "$ROOT" \
 --rna /scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated \
 --model /scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads/Human-GEM-v1.4.1.xml \
 --output data/processed/revised_binary_native_20260910 --native-indicators
