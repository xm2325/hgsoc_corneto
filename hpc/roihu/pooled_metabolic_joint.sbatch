#!/usr/bin/env bash
#SBATCH --job-name=pooled-metabolic-joint
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=384G

set -euo pipefail
: "${SAMPLE_POLICY:?set SAMPLE_POLICY}"
: "${EXPECTED_SAMPLES:?set EXPECTED_SAMPLES}"
: "${OUTPUT_PATH:?set OUTPUT_PATH}"

ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
POOL=${ROOT}/data/processed/rna/pooled_primary_60/tpm
GEM=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads/Human-GEM-v1.4.1.xml
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
cd "${ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH=src
export HGSOC_CORNETO_REPO_COMMIT="$(git rev-parse HEAD)"
export GRB_LICENSE_FILE="${GRB_LICENSE_FILE:-${ROOT}/gurobi_xm.lic}"
test -r "${GRB_LICENSE_FILE}"
test -s data/processed/corneto/pooled_primary_60/metabolic_input_gate.json
test ! -e "${OUTPUT_PATH}"
srun "${PY}" scripts/run_corneto_pooled_metabolic.py \
    --expression "${POOL}/gene_tpm.tsv.gz" \
    --sample-manifest "${POOL}/sample_manifest.tsv" \
    --pooled-receipt "${POOL}/pooled_expression_receipt.json" \
    --human-gem "${GEM}" \
    --output "${OUTPUT_PATH}" \
    --sample-policy "${SAMPLE_POLICY}" \
    --expected-samples "${EXPECTED_SAMPLES}" \
    --max-candidates 25 \
    --growth-fraction 0.9 \
    --joint-lambda 1.0 \
    --solver gurobi
