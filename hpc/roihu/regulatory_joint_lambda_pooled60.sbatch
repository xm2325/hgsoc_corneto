#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=reg-joint-p60

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
REG_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated

LAMBDA_VALUES=(0 0.001 0.01 0.05 0.1 0.25 0.5 1 2)
LAMBDA_LABELS=(l0 l0p001 l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2)
TASK=${SLURM_ARRAY_TASK_ID:?}
LAMBDA=${LAMBDA_VALUES[$TASK]}
LABEL=${LAMBDA_LABELS[$TASK]}
OUT=data/processed/corneto/regulatory_joint_lambda/pooled60/${LABEL}.json
GRAPH=data/processed/corneto/regulatory_joint_lambda/frozen_candidate_graph60.tsv
mkdir -p "$(dirname "$OUT")"
test ! -e "$OUT"

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
  --frozen-graph "$GRAPH" \
  --output "$OUT" --analysis-name pooled60 --primary-only \
  --lambda-reg-reported "$LAMBDA" --max-seconds 258000 --solver gurobi
