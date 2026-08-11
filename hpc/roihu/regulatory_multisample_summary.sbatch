#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=reg-ms-summary

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
GRID=data/processed/corneto/regulatory_multisample/v1/grid
OUT=data/processed/corneto/regulatory_multisample/v1/pooled_vs_cohort_summary.json
LAMBDA_VALUES=(0 0.001 0.01 0.05 0.1 0.25 0.5 1 2)
LAMBDA_LABELS=(l0 l0p001 l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2)
STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)
ARGS=()

for INDEX in "${!LAMBDA_VALUES[@]}"; do
  VALUE=${LAMBDA_VALUES[$INDEX]}
  LABEL=${LAMBDA_LABELS[$INDEX]}
  ARGS+=(--receipt "pooled||${VALUE}|${GRID}/pooled/${LABEL}.json")
  for STUDY in "${STUDIES[@]}"; do
    ARGS+=(--receipt "cohort|${STUDY}|${VALUE}|${GRID}/${STUDY}/${LABEL}.json")
  done
done

test ! -e "$OUT"
exec srun "$PY" scripts/summarize_regulatory_multisample_grid.py \
  "${ARGS[@]}" --output "$OUT"
