#!/usr/bin/env bash
# Response-blind signed-regulon/PKN pilot. CORNETO regulatory solve is explicit
# blocked until a pinned adapter exists; no solver or license is used.
#SBATCH --job-name=hgsoc-regulatory-pilot
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/regulatory-pilot-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna
VENV_ROOT=${HGSOC_CORNETO_ENV:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env}
REGULATORY_SOURCE_ROOT=${REGULATORY_SOURCE_ROOT:-${REPO_ROOT}/data/processed/regulatory_sources}
OUTPUT_ROOT=${REGULATORY_OUTPUT_ROOT:-${REPO_ROOT}/data/processed/regulatory/E-MTAB-14568/pilot_8}
MAX_SAMPLES=${MAX_SAMPLES:-8}

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12
export PYTHONNOUSERSITE=1
export PYTHONPATH=${REPO_ROOT}/src
cd "${REPO_ROOT}"

test -s "${RNA_ROOT}/aggregated/E-MTAB-14568/gene_tpm.tsv.gz"
test -s "${REGULATORY_SOURCE_ROOT}/collectri.tsv"
test -s "${REGULATORY_SOURCE_ROOT}/omnipath.tsv"
test ! -e "${OUTPUT_ROOT}/regulatory_pilot_summary.json"

exec "${VENV_ROOT}/bin/python" scripts/run_corneto_regulatory_pilot.py \
    --expression "${RNA_ROOT}/aggregated/E-MTAB-14568/gene_tpm.tsv.gz" \
    --manifest data/processed/metadata/ocm_master_manifest.tsv \
    --collectri "${REGULATORY_SOURCE_ROOT}/collectri.tsv" \
    --omnipath "${REGULATORY_SOURCE_ROOT}/omnipath.tsv" \
    --output-dir "${OUTPUT_ROOT}" \
    --study E-MTAB-14568 \
    --max-samples "${MAX_SAMPLES}" \
    --top-tfs "${TOP_TFS:-50}" \
    --expression-transform log1p_tpm
