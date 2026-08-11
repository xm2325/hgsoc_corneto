** WARNING: connection is not using a post-quantum key exchange algorithm.
** This session may be vulnerable to "store now, decrypt later" attacks.
** The server may need to be upgraded. See https://openssh.com/pq.html
#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=hgsoc-regulatory

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
STUDY_ACCESSION=${STUDY_ACCESSION:?set STUDY_ACCESSION}
MAX_SAMPLES=${MAX_SAMPLES:?set MAX_SAMPLES}
EXPRESSION=${EXPRESSION:?set absolute expression matrix path}
OUTPUT=${OUTPUT:?set output receipt path relative to repo}
SOURCE_DIR=${SOURCE_DIR:-/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-$ROOT/gurobi_xm.lic}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
test -x "$PY"
test -s "$EXPRESSION"
test -s "$SOURCE_DIR/collectri_tf_target.tsv.gz"
test -s "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz"
test ! -e "$OUTPUT"
exec srun "$PY" scripts/run_corneto_regulatory_pilot.py \
  --expression "$EXPRESSION" \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --collectri "$SOURCE_DIR/collectri_tf_target.tsv.gz" \
  --pkn "$SOURCE_DIR/omnipath_signed_pkn.tsv.gz" \
  --study "$STUDY_ACCESSION" --output "$OUTPUT" --primary-only \
  --max-samples "$MAX_SAMPLES" --max-inputs "${MAX_INPUTS:-3}" \
  --max-outputs "${MAX_OUTPUTS:-6}" --max-depth "${MAX_DEPTH:-3}" \
  --max-edges "${MAX_EDGES:-3000}" --max-flow "${MAX_FLOW:-100}" \
  --lambda-reg "${LAMBDA_REG:-0}" --max-seconds "${MAX_SECONDS:-300}" \
  --solver "${SOLVER:-gurobi}"
