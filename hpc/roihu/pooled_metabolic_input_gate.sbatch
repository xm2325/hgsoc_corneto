#!/usr/bin/env bash
#SBATCH --job-name=pooled-metabolic-gate
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/pooled-metabolic-gate-%j.out

set -euo pipefail
ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
POOL=${ROOT}/data/processed/rna/pooled_primary_60/tpm
GEM=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads/Human-GEM-v1.4.1.xml
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
OUT=${ROOT}/data/processed/corneto/pooled_primary_60/metabolic_input_gate.json
cd "${ROOT}"
test ! -e "${OUT}"
srun "${PY}" scripts/validate_pooled_metabolic_input.py \
    --matrix "${POOL}/gene_tpm.tsv.gz" \
    --sample-manifest "${POOL}/sample_manifest.tsv" \
    --pooled-receipt "${POOL}/pooled_expression_receipt.json" \
    --human-gem "${GEM}" \
    --output "${OUT}"
