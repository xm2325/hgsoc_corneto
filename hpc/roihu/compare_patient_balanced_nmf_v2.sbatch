#!/usr/bin/env bash
# Retry with support for the NMF runner's technical_state assignment column.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=nmf-patient-compare-v2

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
POOLED=data/processed/rna/pooled_primary60/nmf/rank_3/assignments.tsv
BALANCED=data/processed/rna/pooled_patient_balanced52/nmf/rank_3/assignments.tsv
OUT=data/processed/rna/pooled_patient_balanced52/patient_balanced_vs_pooled_rank3_v2.json
test -s "$POOLED"; test -s "$BALANCED"; test ! -e "$OUT"
exec srun "$PY" scripts/compare_patient_balanced_nmf.py \
  --pooled "$POOLED" --patient-balanced "$BALANCED" --output "$OUT"
