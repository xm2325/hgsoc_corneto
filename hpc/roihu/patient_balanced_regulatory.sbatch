#!/usr/bin/env bash
# Patient-balanced true-joint regulatory sensitivity at lambda=0.001.
# One solver task only; it does not wait for metabolic baselines.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --job-name=reg-patient-balanced

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}

BASE=data/processed/corneto/regulatory_multisample/patient_balanced52
BUNDLE="$BASE/input_bundle.json"
SELECTED=data/processed/rna/pooled_patient_balanced52/input/selected_samples.tsv
OUT="$BASE/regulatory_l0p001.json"
test ! -e "$BASE"
test -s data/processed/corneto/regulatory_multisample/v1/input_bundle.json
test -s "$SELECTED"
mkdir -p "$BASE"
srun "$PY" scripts/prepare_patient_balanced_regulatory_bundle.py \
  --bundle data/processed/corneto/regulatory_multisample/v1/input_bundle.json \
  --selected "$SELECTED" --output "$BUNDLE"
exec srun "$PY" scripts/run_corneto_regulatory_multisample.py \
  --bundle "$BUNDLE" --output "$OUT" --mode pooled \
  --lambda-reg 0.001 --lambda-scaling mean_fit --solver gurobi \
  --max-flow 100 --max-seconds 258000
