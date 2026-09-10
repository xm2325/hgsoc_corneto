#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --job-name=revised-gpr-pfba
set -euo pipefail
ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
cd "$ROOT"
export PYTHONPATH="$ROOT/scripts:$ROOT/src" PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
exec srun /scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python \
 scripts/validate_revised_metabolic_model.py --root "$ROOT" \
 --rna-root /scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated \
 --model /scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads/Human-GEM-v1.4.1.xml \
 --output data/processed/revised_gpr_pfba_20260910
