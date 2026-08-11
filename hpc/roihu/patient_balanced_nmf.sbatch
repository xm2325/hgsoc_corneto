#!/usr/bin/env bash
# Patient-balanced (one OCM per patient) NMF sensitivity.
# This is independent of the metabolic CORNETO jobs.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=nmf-patient-balanced

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src" OMP_NUM_THREADS=8 MKL_NUM_THREADS=8

BASE=data/processed/rna/pooled_patient_balanced52
INPUT="$BASE/input"
OUT="$BASE/nmf"
MATRIX=data/processed/rna/pooled_primary_60/log1p_tpm/gene_log1p_tpm.tsv.gz
MANIFEST=data/processed/metadata/ocm_master_manifest.tsv

test -s "$MATRIX"
test -s "$MANIFEST"
test ! -e "$INPUT"
test ! -e "$OUT"

srun "$PY" scripts/prepare_patient_balanced_nmf.py \
  --matrix "$MATRIX" --manifest "$MANIFEST" --output-dir "$INPUT"

exec srun "$PY" scripts/run_nmf_primary.py \
  --analysis-id pooled-patient-balanced-52 \
  --cohort "patient_balanced_primary::$INPUT/gene_log1p_tpm_patient_balanced.tsv.gz::$INPUT/sample_qc.tsv" \
  --output-dir "$OUT" --top-genes 6000 --ranks 2,3,4,5,6,7,8,9,10 \
  --starts 100 --seed-base 605200
