#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --job-name=reg-ms-grid

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
BUNDLE=data/processed/corneto/regulatory_multisample/v1/input_bundle.json

LAMBDA_VALUES=(0 0.001 0.01 0.05 0.1 0.25 0.5 1 2)
LAMBDA_LABELS=(l0 l0p001 l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2)
STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)
TASK=${SLURM_ARRAY_TASK_ID:?}
GROUP_INDEX=$((TASK / 9))
LAMBDA_INDEX=$((TASK % 9))
LAMBDA=${LAMBDA_VALUES[$LAMBDA_INDEX]}
LABEL=${LAMBDA_LABELS[$LAMBDA_INDEX]}
ARGS=(--mode pooled)
GROUP=pooled
if (( GROUP_INDEX > 0 )); then
  STUDY=${STUDIES[$((GROUP_INDEX - 1))]}
  GROUP=$STUDY
  ARGS=(--mode cohort --study "$STUDY")
fi
OUT=data/processed/corneto/regulatory_multisample/v1/grid/${GROUP}/${LABEL}.json

test -s "$BUNDLE"
test ! -e "$OUT"
exec srun "$PY" scripts/run_corneto_regulatory_multisample.py \
  --bundle "$BUNDLE" --output "$OUT" "${ARGS[@]}" \
  --lambda-reg "$LAMBDA" --lambda-scaling mean_fit --solver gurobi \
  --max-flow 100 --max-seconds 258000
