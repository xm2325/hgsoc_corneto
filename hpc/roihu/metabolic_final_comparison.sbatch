#!/usr/bin/env bash
#SBATCH --job-name=metabolic-final-compare
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/metabolic-final-compare-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
PYTHON_BIN=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
OUTPUT_PATH="${OUTPUT_PATH:-data/processed/corneto/metabolic_full_cohort_comparison_v2.json}"

cd "${REPO_ROOT}"
test ! -e "${OUTPUT_PATH}"
srun "${PYTHON_BIN}" scripts/validate_corneto_cohort_receipts.py \
    --receipt E-MTAB-7223=data/processed/corneto/E-MTAB-7223/full_direct_b25.json \
    --receipt E-MTAB-10801=data/processed/corneto/E-MTAB-10801/full_direct_b25.json \
    --receipt E-MTAB-11000=data/processed/corneto/E-MTAB-11000/full_direct_b25.json \
    --receipt E-MTAB-14568=data/processed/corneto/E-MTAB-14568/full_direct_b25.json \
    --output "${OUTPUT_PATH}"
