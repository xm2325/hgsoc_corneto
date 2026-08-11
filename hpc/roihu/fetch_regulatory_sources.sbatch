#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=hgsoc-reg-sources

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
OUT=${REGULATORY_OUT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
test -x "$PY"
test ! -e "$OUT/regulatory_source_receipt.json"
exec "$PY" scripts/fetch_regulatory_sources.py --out-dir "$OUT"
