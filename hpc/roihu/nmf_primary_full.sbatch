#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=nmf-primary-full

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
RNA_ROOT=${RNA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src" OMP_NUM_THREADS=8 MKL_NUM_THREADS=8

case "${SLURM_ARRAY_TASK_ID:?}" in
  0) TARGET=E-MTAB-7223; STUDIES=(E-MTAB-7223); SEED=722300 ;;
  1) TARGET=E-MTAB-10801; STUDIES=(E-MTAB-10801); SEED=1080100 ;;
  2) TARGET=E-MTAB-11000; STUDIES=(E-MTAB-11000); SEED=1100000 ;;
  3) TARGET=E-MTAB-14568; STUDIES=(E-MTAB-14568); SEED=1456800 ;;
  4) TARGET=pooled60; STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568); SEED=6000000 ;;
  *) echo "Unknown array index: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac

OUT="data/processed/rna/nmf_primary_v2/${TARGET}"
if [[ -e "$OUT" ]]; then
  test -s "$OUT/nmf_receipt.json"
  exit 0
fi
COHORT_ARGS=()
for STUDY in "${STUDIES[@]}"; do
  MATRIX="$RNA_ROOT/$STUDY/gene_log1p_tpm.tsv.gz"
  QC="data/processed/rna/$STUDY/aggregation/sample_qc.tsv"
  test -s "$MATRIX"
  test -s "$QC"
  COHORT_ARGS+=(--cohort "$STUDY::$MATRIX::$QC")
done
exec srun "$PY" scripts/run_nmf_primary.py \
  --analysis-id "primary-${TARGET}" "${COHORT_ARGS[@]}" --output-dir "$OUT" \
  --top-genes 6000 --ranks 2,3,4,5,6,7,8,9,10 --starts 100 --seed-base "$SEED"
