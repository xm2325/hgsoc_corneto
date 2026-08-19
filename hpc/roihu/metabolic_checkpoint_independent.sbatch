#!/bin/bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

: "${STUDY_ACCESSION:?set STUDY_ACCESSION}"
: "${SLURM_ARRAY_TASK_ID:?array task required}"

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
INPUT_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads
PYTHON=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
CHECKPOINT_ROOT="${REPO_ROOT}/data/processed/corneto/${STUDY_ACCESSION}/checkpoint_b25"

cd "${REPO_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH=src:scripts
mkdir -p "${CHECKPOINT_ROOT}/independent"

srun "${PYTHON}" scripts/run_corneto_metabolic_checkpoint.py independent \
    --context "${CHECKPOINT_ROOT}/context.json" \
    --human-gem "${INPUT_ROOT}/Human-GEM-v1.4.1.xml" \
    --array-index "${SLURM_ARRAY_TASK_ID}" \
    --output-dir "${CHECKPOINT_ROOT}/independent"
