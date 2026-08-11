#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=pooled60-gate

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
RNA_ROOT=${RNA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python

case "${SLURM_ARRAY_TASK_ID:?}" in
  0) VALUE=tpm ;;
  1) VALUE=log1p_tpm ;;
  *) echo "unknown array index: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac

OUT=data/processed/rna/pooled_primary_60/${VALUE}
test ! -e "$OUT"
exec srun "$PY" scripts/build_pooled_primary_expression.py \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --matrix E-MTAB-7223="$RNA_ROOT/E-MTAB-7223/gene_${VALUE}.tsv.gz" \
  --matrix E-MTAB-10801="$RNA_ROOT/E-MTAB-10801/gene_${VALUE}.tsv.gz" \
  --matrix E-MTAB-11000="$RNA_ROOT/E-MTAB-11000/gene_${VALUE}.tsv.gz" \
  --matrix E-MTAB-14568="$RNA_ROOT/E-MTAB-14568/gene_${VALUE}.tsv.gz" \
  --output-dir "$OUT" --value-name "$VALUE" \
  --expected-samples 60 --expected-patients 52
