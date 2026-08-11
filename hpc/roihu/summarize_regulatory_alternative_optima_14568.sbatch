#!/bin/bash
#SBATCH --job-name=reg-alt-summary
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/reg-alt-summary-%j.out

set -euo pipefail

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
INPUT=${ROOT}/data/processed/corneto/alternative_optima/E-MTAB-14568/l0p1
RETRY_INPUT=${RETRY_INPUT:-}
OUTPUT=${OUTPUT:-${ROOT}/data/processed/corneto/alternative_optima/E-MTAB-14568/l0p1_summary.json}

cd "$ROOT"
test -x "$PY"
test ! -e "$OUTPUT"

INPUT_ARGS=(--input-dir "$INPUT")
if [[ -n "$RETRY_INPUT" ]]; then
  INPUT_ARGS+=(--input-dir "$RETRY_INPUT")
fi

exec srun "$PY" scripts/summarize_regulatory_alternative_optima.py \
  "${INPUT_ARGS[@]}" \
  --expected-samples 27 \
  --study E-MTAB-14568 \
  --output "$OUTPUT"
