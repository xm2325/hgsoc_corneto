#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=reg-joint-graph

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
REG_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated
OUT=data/processed/corneto/regulatory_joint_lambda/frozen_candidate_graph60_receipt.json
GRAPH=data/processed/corneto/regulatory_joint_lambda/frozen_candidate_graph60.tsv
mkdir -p "$(dirname "$OUT")"
test ! -e "$OUT"
test ! -e "$GRAPH"

exec srun "$PY" scripts/run_corneto_regulatory_joint.py \
  --expression E-MTAB-7223="$RNA_ROOT/E-MTAB-7223/gene_tpm.tsv.gz" \
  --expression E-MTAB-10801="$RNA_ROOT/E-MTAB-10801/gene_tpm.tsv.gz" \
  --expression E-MTAB-11000="$RNA_ROOT/E-MTAB-11000/gene_tpm.tsv.gz" \
  --expression E-MTAB-14568="$RNA_ROOT/E-MTAB-14568/gene_tpm.tsv.gz" \
  --expected-count E-MTAB-7223=9 --expected-count E-MTAB-10801=13 \
  --expected-count E-MTAB-11000=11 --expected-count E-MTAB-14568=27 \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --collectri "$REG_ROOT/collectri_tf_target.tsv.gz" \
  --pkn "$REG_ROOT/omnipath_signed_pkn.tsv.gz" \
  --frozen-graph "$GRAPH" --prepare-only \
  --output "$OUT" --analysis-name frozen_graph_all60 --primary-only
