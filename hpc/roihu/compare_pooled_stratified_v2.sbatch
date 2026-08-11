#!/usr/bin/env bash
# Final pooled-vs-cohort post-processing for the normalized true-joint grid.
# This job reads completed receipts only; it never invokes a solver.
# The v1 wrapper pointed at the legacy regulatory_joint_lambda tree. The
# normalized 591595 grid lives under regulatory_multisample/v1/grid instead.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=pooled-stratified-v2

set -euo pipefail

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"

GRID=data/processed/corneto/regulatory_multisample/v1/grid
MANIFEST=data/processed/metadata/ocm_master_manifest.tsv
OUT=data/processed/corneto/pooled_vs_cohort_joint_nmf_comparison_v2.json
TABLES=data/processed/corneto/pooled_vs_cohort_joint_nmf_tables_v2

test -s "$MANIFEST"
test ! -e "$OUT"

LAMBDA_VALUES=(0 0.001 0.01 0.05 0.1 0.25 0.5 1 2)
LAMBDA_LABELS=(l0 l0p001 l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2)
STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)

ARGS=(
  --manifest "$MANIFEST"
  --require-joint-pooled
  --pooled-nmf data/processed/rna/pooled_primary60/nmf/rank_3/assignments.tsv
  --output-json "$OUT"
  --output-dir "$TABLES"
)
test -s "data/processed/rna/pooled_primary60/nmf/rank_3/assignments.tsv"

for index in "${!LAMBDA_VALUES[@]}"; do
  value=${LAMBDA_VALUES[$index]}
  label=${LAMBDA_LABELS[$index]}
  pooled="$GRID/pooled/$label.json"
  test -s "$pooled"
  ARGS+=(--pooled-regulatory "${value}=${pooled}")
  for study in "${STUDIES[@]}"; do
    cohort="$GRID/$study/$label.json"
    test -s "$cohort"
    ARGS+=(--cohort-regulatory "${study}|${value}|${cohort}")
  done
done

for study in "${STUDIES[@]}"; do
  assignments="data/processed/rna/$study/nmf_primary/rank_3/assignments.tsv"
  test -s "$assignments"
  ARGS+=(--cohort-nmf "${study}|${assignments}")
done

exec srun "$PY" scripts/compare_pooled_stratified.py "${ARGS[@]}"
