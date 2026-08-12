#!/usr/bin/env bash
#SBATCH --job-name=metabolic-full-retry
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail

: "${STUDY_ACCESSION:?set STUDY_ACCESSION}"
: "${EXPECTED_SAMPLES:?set EXPECTED_SAMPLES}"

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated
INPUT_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads
PYTHON_BIN=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
OUTPUT_PATH="${OUTPUT_PATH:-${REPO_ROOT}/data/processed/corneto/${STUDY_ACCESSION}/full_direct_b25.json}"
EXPRESSION_PATH="${RNA_ROOT}/${STUDY_ACCESSION}/gene_tpm.tsv.gz"
HUMAN_GEM_PATH="${INPUT_ROOT}/Human-GEM-v1.4.1.xml"

cd "${REPO_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH=src
export HGSOC_CORNETO_REPO_COMMIT="$(git rev-parse HEAD)"
export GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-${REPO_ROOT}/gurobi_xm.lic}"

test -r "${GRB_LICENSE_FILE}"
test -s "${EXPRESSION_PATH}"
test -s "${HUMAN_GEM_PATH}"

# A speculative retry is safe to queue before the original job terminates.
# If the original produced a valid receipt, validate it and exit without
# overwriting.  An existing invalid receipt is retained and fails closed.
if [[ -e "${OUTPUT_PATH}" ]]; then
    "${PYTHON_BIN}" scripts/validate_corneto_14568_pilot.py \
        "${OUTPUT_PATH}" \
        --expected-solver gurobi \
        --expected-samples "${EXPECTED_SAMPLES}" \
        --expected-candidates 25 \
        --expected-growth-fraction 0.9
    printf 'retry_status=skipped_valid_existing_receipt\n'
    exit 0
fi

srun "${PYTHON_BIN}" scripts/run_corneto_14568_pilot.py \
    --study "${STUDY_ACCESSION}" \
    --expression "${EXPRESSION_PATH}" \
    --manifest data/processed/metadata/ocm_master_manifest.tsv \
    --human-gem "${HUMAN_GEM_PATH}" \
    --output "${OUTPUT_PATH}" \
    --primary-only \
    --max-samples "${EXPECTED_SAMPLES}" \
    --max-candidates 25 \
    --expression-transform log1p_tpm \
    --solver gurobi \
    --growth-fraction 0.9 \
    --independent-lambda 0.1 \
    --joint-lambda 1.0
