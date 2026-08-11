#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=nmf-cross-cohort

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python

case "${SLURM_ARRAY_TASK_ID:?}" in
  0) STUDY=E-MTAB-7223 ;;
  1) STUDY=E-MTAB-10801 ;;
  2) STUDY=E-MTAB-11000 ;;
  *) echo "unknown array index: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac
MATRIX=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated/${STUDY}/gene_log1p_tpm.tsv.gz
QC=data/processed/rna/${STUDY}/aggregation/sample_qc.tsv
OUT=data/processed/rna/${STUDY}/nmf_cross_cohort
test -s "$MATRIX"
test -s "$QC"
test ! -e "$OUT"
exec srun "$PY" scripts/run_nmf_cross_cohort.py \
  --study "$STUDY" --matrix "$MATRIX" --sample-qc "$QC" --output-dir "$OUT" \
  --top-genes 6000 --starts 100 --max-rank 8
