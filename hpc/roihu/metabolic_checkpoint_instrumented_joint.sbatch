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

ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
GEM=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads/Human-GEM-v1.4.1.xml
CHECKPOINT_ROOT="${ROOT}/data/processed/corneto/${STUDY_ACCESSION}/checkpoint_b25"
MAX_SECONDS=${MAX_SECONDS:-252000}
MIP_GAP=${MIP_GAP:-0.0001}

cd "${ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH=src:scripts
export GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-${ROOT}/gurobi_xm.lic}"
test -r "${GRB_LICENSE_FILE}"

srun "${PY}" scripts/run_corneto_metabolic_instrumented.py joint \
    --context "${CHECKPOINT_ROOT}/context.json" \
    --human-gem "${GEM}" \
    --output "${CHECKPOINT_ROOT}/joint.json" \
    --attempt-dir "${CHECKPOINT_ROOT}/instrumented_attempts" \
    --artifact-dir "${CHECKPOINT_ROOT}/solver_artifacts" \
    --max-seconds "${MAX_SECONDS}" \
    --mip-gap "${MIP_GAP}" \
    --threads "${SLURM_CPUS_PER_TASK}"
