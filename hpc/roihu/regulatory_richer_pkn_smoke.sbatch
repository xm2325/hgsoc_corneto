#!/usr/bin/env bash
# Two-condition solver smoke for the richer PKN bundle.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=reg-richer-smoke

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}

BUNDLE=data/processed/corneto/regulatory_multisample/richer_v1/input_bundle.json
OUT=data/processed/corneto/regulatory_multisample/richer_v1/smoke/two_condition_l0p001.json
test -s "$BUNDLE"
test ! -e "$OUT"
exec srun "$PY" scripts/run_corneto_regulatory_multisample.py \
  --bundle "$BUNDLE" --output "$OUT" --mode pooled --max-conditions 2 \
  --lambda-reg 0.001 --lambda-scaling mean_fit --solver gurobi \
  --max-flow 100 --max-seconds 600
