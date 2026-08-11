#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=reg-lambda-audit

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python

STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)
COUNTS=(9 13 11 27)
VALUES=(0.01 0.05 0.1 0.25 0.5 1 2)
LABELS=(l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2)
INCLUDE_L0P001=${INCLUDE_L0P001:-0}
ARGS=()

for index in "${!STUDIES[@]}"; do
  study=${STUDIES[$index]}
  count=${COUNTS[$index]}
  baseline="data/processed/corneto/${study}/regulatory_full_gurobi_${count}_l0.json"
  test -s "$baseline"
  ARGS+=(--receipt "${study}|0|${baseline}")
  if [[ "$INCLUDE_L0P001" == 1 ]]; then
    path="data/processed/corneto/lambda_grid_all60/${study}/l0p001.json"
    test -s "$path"
    ARGS+=(--receipt "${study}|0.001|${path}")
  fi
  for lambda_index in "${!VALUES[@]}"; do
    value=${VALUES[$lambda_index]}
    label=${LABELS[$lambda_index]}
    path="data/processed/corneto/lambda_grid_all60/${study}/${label}.json"
    test -s "$path"
    ARGS+=(--receipt "${study}|${value}|${path}")
  done
done

if [[ "$INCLUDE_L0P001" == 1 ]]; then
  OUT=data/processed/corneto/regulatory_lambda_all60_audit_grid9.json
else
  OUT=data/processed/corneto/regulatory_lambda_all60_audit_grid8.json
fi
test ! -e "$OUT"
exec srun "$PY" scripts/summarize_regulatory_lambda_all60.py "${ARGS[@]}" --output "$OUT"
