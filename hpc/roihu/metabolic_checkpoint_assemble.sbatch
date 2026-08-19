#!/bin/bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

: "${STUDY_ACCESSION:?set STUDY_ACCESSION}"

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
PYTHON=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
CHECKPOINT_ROOT="${REPO_ROOT}/data/processed/corneto/${STUDY_ACCESSION}/checkpoint_b25"
FINAL_OUTPUT="${REPO_ROOT}/data/processed/corneto/${STUDY_ACCESSION}/full_direct_b25.json"

cd "${REPO_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH=src:scripts

srun "${PYTHON}" scripts/run_corneto_metabolic_checkpoint.py assemble \
    --context "${CHECKPOINT_ROOT}/context.json" \
    --independent-dir "${CHECKPOINT_ROOT}/independent" \
    --joint-receipt "${CHECKPOINT_ROOT}/joint.json" \
    --output "${FINAL_OUTPUT}"
