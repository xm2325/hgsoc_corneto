#!/usr/bin/env bash
# One-at-a-time richer-PKN λ=0.001 sensitivity: pooled + four cohorts.
# The %1 throttle deliberately leaves Gurobi session and node memory headroom.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --job-name=reg-richer-grid

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}

BUNDLE=data/processed/corneto/regulatory_multisample/richer_v1/input_bundle.json
test -s "$BUNDLE"
TASK=${SLURM_ARRAY_TASK_ID:?}
GROUP=pooled
ARGS=(--mode pooled)
if (( TASK > 0 )); then
  STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)
  STUDY=${STUDIES[$((TASK - 1))]}
  GROUP=$STUDY
  ARGS=(--mode cohort --study "$STUDY")
fi
OUT=data/processed/corneto/regulatory_multisample/richer_v1/grid/$GROUP/l0p001.json
test ! -e "$OUT"
exec srun "$PY" scripts/run_corneto_regulatory_multisample.py \
  --bundle "$BUNDLE" --output "$OUT" "${ARGS[@]}" \
  --lambda-reg 0.001 --lambda-scaling mean_fit --solver gurobi \
  --max-flow 100 --max-seconds 258000
