#!/usr/bin/env bash
# Response-blind longitudinal edge-change summary across all four receipts.
# It compares primary family runs only and does not use phenotype.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=reg-longitudinal-joint

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"

OUT=data/processed/corneto/regulatory_longitudinal_joint_l0p001.json
FAMILIES=data/processed/taylor_preflight/acquired_family_design.tsv
test ! -e "$OUT"
test -s "$FAMILIES"

ARGS=(--families "$FAMILIES" --output "$OUT" --lambda-nominal 0.001)
for study in E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568; do
  receipt=data/processed/corneto/regulatory_multisample/v1/grid/$study/l0p001.json
  test -s "$receipt"
  ARGS+=(--receipt "$receipt")
done
exec srun "$PY" scripts/analyze_regulatory_longitudinal_joint.py "${ARGS[@]}"
