#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=ext-gse208-input
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-gse208-input-%j.out

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
EXT=${EXTERNAL_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}
COUNTS=${EXT}/GSE208216/GSE208216_ViasMorrill_RNASeq_raw_counts.txt.gz
GENE_MAP=${RNA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated}/E-MTAB-14568/gene_metadata.tsv
OUT=${EXT}/GSE208216/corneto_input_v1

cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME
test -s "${COUNTS}"
test -s "${GENE_MAP}"
test ! -e "${OUT}"
srun "${PY}" scripts/prepare_gse208216_corneto_input.py \
    --counts "${COUNTS}" --gene-map "${GENE_MAP}" \
    --contract config/external_validation/gse208216.json --output-dir "${OUT}"

"${PY}" -c '
import json, sys
r=json.load(open(sys.argv[1]))
assert r["status"] == "completed"
assert r["dimensions"]["samples"] == 14
assert r["dimensions"]["pdo_models"] == 11
assert r["dimensions"]["fallopian_tube_models"] == 3
' "${OUT}/receipt.json"
