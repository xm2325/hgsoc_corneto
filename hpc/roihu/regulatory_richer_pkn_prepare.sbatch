#!/usr/bin/env bash
# Build-only richer bounded PKN bundle; no solver/license checkout.
# Compared with v1: more inputs/outputs and one extra upstream BFS level.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --job-name=reg-richer-input

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
RNA_ROOT=${RNA_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated}
SOURCE_DIR=${REGULATORY_SOURCE_DIR:-/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"

OUT=data/processed/corneto/regulatory_multisample/richer_v1/input_bundle.json
test ! -e "$OUT"
exec srun "$PY" scripts/prepare_regulatory_multisample.py \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --expression E-MTAB-7223="$RNA_ROOT/E-MTAB-7223/gene_tpm.tsv.gz" \
  --expression E-MTAB-10801="$RNA_ROOT/E-MTAB-10801/gene_tpm.tsv.gz" \
  --expression E-MTAB-11000="$RNA_ROOT/E-MTAB-11000/gene_tpm.tsv.gz" \
  --expression E-MTAB-14568="$RNA_ROOT/E-MTAB-14568/gene_tpm.tsv.gz" \
  --collectri "$SOURCE_DIR/collectri_tf_target.tsv.gz" \
  --pkn "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz" \
  --output "$OUT" --min-targets 5 --max-inputs 5 --max-outputs 10 \
  --max-depth 4 --max-edges 10000
