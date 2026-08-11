#!/usr/bin/env bash
# Solver-free comparison after patient-balanced regulatory solve.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=reg-patient-compare

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
POOLED=data/processed/corneto/regulatory_multisample/v1/grid/pooled/l0p001.json
BALANCED=data/processed/corneto/regulatory_multisample/patient_balanced52/regulatory_l0p001.json
OUT=data/processed/corneto/regulatory_multisample/patient_balanced52/pooled_vs_patient_balanced.json
test -s "$POOLED"; test -s "$BALANCED"; test ! -e "$OUT"
exec srun "$PY" scripts/compare_patient_balanced_regulatory.py \
  --pooled "$POOLED" --patient-balanced "$BALANCED" --output "$OUT"
