#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=reg-lambda-all60

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}

STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)
SAMPLE_COUNTS=(9 13 11 27)
LAMBDA_VALUES=(0.01 0.05 0.1 0.25 0.5 1 2)
LAMBDA_LABELS=(l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2)
TASK=${SLURM_ARRAY_TASK_ID:?}
COHORT_INDEX=$((TASK / 7))
LAMBDA_INDEX=$((TASK % 7))
STUDY=${STUDIES[$COHORT_INDEX]}
N=${SAMPLE_COUNTS[$COHORT_INDEX]}
LAMBDA=${LAMBDA_VALUES[$LAMBDA_INDEX]}
LABEL=${LAMBDA_LABELS[$LAMBDA_INDEX]}
EXPRESSION=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated/${STUDY}/gene_tpm.tsv.gz
SOURCE_DIR=/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory
OUT=data/processed/corneto/lambda_grid_all60/${STUDY}/${LABEL}.json
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
mkdir -p "$(dirname "$OUT")"
test ! -e "$OUT"
test -s "$EXPRESSION"
test -s "$SOURCE_DIR/collectri_tf_target.tsv.gz"
test -s "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz"
exec srun "$PY" scripts/run_corneto_regulatory_pilot.py \
  --expression "$EXPRESSION" \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --collectri "$SOURCE_DIR/collectri_tf_target.tsv.gz" \
  --pkn "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz" \
  --study "$STUDY" --output "$OUT" --primary-only \
  --max-samples "$N" --max-inputs 3 --max-outputs 6 \
  --max-depth 3 --max-edges 3000 --max-flow 100 \
  --lambda-reg "$LAMBDA" --max-seconds 300 --solver gurobi
