#!/bin/bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

: "${STUDY_ACCESSION:?set STUDY_ACCESSION}"
: "${EXPECTED_SAMPLES:?set EXPECTED_SAMPLES}"

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated
INPUT_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads
PYTHON=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
CHECKPOINT_ROOT="${REPO_ROOT}/data/processed/corneto/${STUDY_ACCESSION}/checkpoint_b25"

cd "${REPO_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH=src:scripts
export GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-${REPO_ROOT}/gurobi_xm.lic}"
test -r "${GRB_LICENSE_FILE}"
export HGSOC_CORNETO_REPO_COMMIT
HGSOC_CORNETO_REPO_COMMIT=$(git rev-parse HEAD)
mkdir -p "${CHECKPOINT_ROOT}"

srun "${PYTHON}" scripts/run_corneto_metabolic_checkpoint.py prepare \
    --expression "${RNA_ROOT}/${STUDY_ACCESSION}/gene_tpm.tsv.gz" \
    --manifest data/processed/metadata/ocm_master_manifest.tsv \
    --human-gem "${INPUT_ROOT}/Human-GEM-v1.4.1.xml" \
    --study "${STUDY_ACCESSION}" \
    --expected-samples "${EXPECTED_SAMPLES}" \
    --output "${CHECKPOINT_ROOT}/context.json"
