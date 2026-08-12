#!/usr/bin/env bash
#SBATCH --job-name=tpi1-model-gate
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/tpi1-model-gate-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
MODEL=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads/Human-GEM-v1.4.1.xml
PYTHON_BIN=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
OUTPUT_PATH="${OUTPUT_PATH:-data/processed/corneto/tpi1_model_gate.json}"

cd "${REPO_ROOT}"
test -s "${MODEL}"
test ! -e "${OUTPUT_PATH}"
srun "${PYTHON_BIN}" scripts/check_tpi1_model_gate.py \
    --sbml "${MODEL}" \
    --output "${OUTPUT_PATH}"
