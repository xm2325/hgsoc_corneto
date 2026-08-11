#!/usr/bin/env bash
# Solver-free comparison of narrow vs richer PKN receipts.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=reg-richer-compare

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
OUT=data/processed/corneto/regulatory_multisample/richer_v1/narrow_vs_richer_l0p001.json
test ! -e "$OUT"
for group in pooled E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568; do
  test -s "data/processed/corneto/regulatory_multisample/v1/grid/$group/l0p001.json"
  test -s "data/processed/corneto/regulatory_multisample/richer_v1/grid/$group/l0p001.json"
done
exec srun "$PY" scripts/compare_richer_pkn_sensitivity.py \
  --narrow-root data/processed/corneto/regulatory_multisample/v1/grid \
  --richer-root data/processed/corneto/regulatory_multisample/richer_v1/grid \
  --output "$OUT"
