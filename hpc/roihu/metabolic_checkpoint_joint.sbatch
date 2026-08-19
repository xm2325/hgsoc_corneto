#!/bin/bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=196G

set -euo pipefail

: "${STUDY_ACCESSION:?set STUDY_ACCESSION}"

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
INPUT_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads
PYTHON=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
CHECKPOINT_ROOT="${REPO_ROOT}/data/processed/corneto/${STUDY_ACCESSION}/checkpoint_b25"

cd "${REPO_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH=src:scripts

srun "${PYTHON}" scripts/run_corneto_metabolic_checkpoint.py joint \
    --context "${CHECKPOINT_ROOT}/context.json" \
    --human-gem "${INPUT_ROOT}/Human-GEM-v1.4.1.xml" \
    --output "${CHECKPOINT_ROOT}/joint.json"
