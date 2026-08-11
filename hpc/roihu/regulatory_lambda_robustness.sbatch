#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=reg-lambda-full27

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}

case "${SLURM_ARRAY_TASK_ID:?}" in
  0) LAMBDA=0.01; LABEL=l0p01 ;;
  1) LAMBDA=0.05; LABEL=l0p05 ;;
  2) LAMBDA=0.1;  LABEL=l0p1  ;;
  3) LAMBDA=0.25; LABEL=l0p25 ;;
  4) LAMBDA=0.5;  LABEL=l0p5  ;;
  5) LAMBDA=1;    LABEL=l1    ;;
  6) LAMBDA=2;    LABEL=l2    ;;
  *) echo "unknown array index: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac

OUT="$ROOT/data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_${LABEL}.json"
EXPRESSION=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated/E-MTAB-14568/gene_tpm.tsv.gz
SOURCE_DIR=/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
test ! -e "$OUT"
test -s "$EXPRESSION"
test -s "$SOURCE_DIR/collectri_tf_target.tsv.gz"
test -s "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz"

exec srun "$PY" scripts/run_corneto_regulatory_pilot.py \
  --expression "$EXPRESSION" \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --collectri "$SOURCE_DIR/collectri_tf_target.tsv.gz" \
  --pkn "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz" \
  --study E-MTAB-14568 --output "$OUT" --primary-only \
  --max-samples 27 --max-inputs 3 --max-outputs 6 \
  --max-depth 3 --max-edges 3000 --max-flow 100 \
  --lambda-reg "$LAMBDA" --max-seconds 300 --solver gurobi
