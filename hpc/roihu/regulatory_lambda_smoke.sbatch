#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=reg-lambda-smoke

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}
case "${SLURM_ARRAY_TASK_ID:?}" in
  0) STUDY=E-MTAB-7223;  N=9  ;;
  1) STUDY=E-MTAB-10801; N=13 ;;
  2) STUDY=E-MTAB-11000; N=11 ;;
  3) STUDY=E-MTAB-14568; N=27 ;;
  *) echo "unknown array index: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac
EXPRESSION=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated/${STUDY}/gene_tpm.tsv.gz
SOURCE_DIR=/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory
OUT=data/processed/corneto/lambda_grid_smoke/${STUDY}.json
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
  --max-samples 1 --max-inputs 3 --max-outputs 6 \
  --max-depth 3 --max-edges 3000 --max-flow 100 \
  --lambda-reg 0.25 --max-seconds 300 --solver gurobi
