#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=pooled-stratified-audit

set -euo pipefail

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1

LAMBDA_VALUES=(0 0.001 0.01 0.05 0.1 0.25 0.5 1 2)
LAMBDA_LABELS=(l0 l0p001 l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2)
STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)

ARGS=(
  --manifest data/processed/metadata/ocm_master_manifest.tsv
  --require-joint-pooled
  --pooled-nmf data/processed/rna/pooled_primary60/nmf/rank_3/assignments.tsv
  --output-json data/processed/corneto/pooled_vs_cohort_joint_nmf_comparison.json
  --output-dir data/processed/corneto/pooled_vs_cohort_joint_nmf_tables
)

for index in "${!LAMBDA_VALUES[@]}"; do
  value=${LAMBDA_VALUES[$index]}
  label=${LAMBDA_LABELS[$index]}
  pooled=data/processed/corneto/regulatory_joint_lambda/pooled60/${label}.json
  test -s "$pooled"
  ARGS+=(--pooled-regulatory "${value}=${pooled}")
  for study in "${STUDIES[@]}"; do
    cohort=data/processed/corneto/regulatory_joint_lambda/cohorts/${study}/${label}.json
    test -s "$cohort"
    ARGS+=(--cohort-regulatory "${study}|${value}|${cohort}")
  done
done

for study in "${STUDIES[@]}"; do
  assignments=data/processed/rna/${study}/nmf_primary/rank_3/assignments.tsv
  test -s "$assignments"
  ARGS+=(--cohort-nmf "${study}|${assignments}")
done

srun "$PY" scripts/compare_pooled_stratified.py "${ARGS[@]}"

# Rank 2 was more stable than rank 3 in the completed 100-start screening for
# all five targets.  Compare it separately so the frozen rank-3 benchmark is
# retained without treating k=3 as the only defensible structure.
NMF2_ARGS=(
  --manifest data/processed/metadata/ocm_master_manifest.tsv
  --pooled-nmf data/processed/rna/pooled_primary60/nmf/rank_2/assignments.tsv
  --output-json data/processed/corneto/pooled_vs_cohort_nmf_rank2_comparison.json
  --output-dir data/processed/corneto/pooled_vs_cohort_nmf_rank2_tables
)
for study in "${STUDIES[@]}"; do
  assignments=data/processed/rna/${study}/nmf_primary/rank_2/assignments.tsv
  test -s "$assignments"
  NMF2_ARGS+=(--cohort-nmf "${study}|${assignments}")
done
srun "$PY" scripts/compare_pooled_stratified.py "${NMF2_ARGS[@]}"
