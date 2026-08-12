#!/usr/bin/env bash
# Fail-closed, non-Gurobi smoke for a prepared external expression bundle.
# Required exports: STUDY_ACCESSION, EXPRESSION, MANIFEST, OUTPUT.
# This is a computational/input gate, not an external biological result.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=ext-reg-smoke
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-reg-smoke-%j.out

set -euo pipefail
umask 027

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
SOURCE_DIR=${SOURCE_DIR:-/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory}
: "${STUDY_ACCESSION:?Set external study accession}"
: "${EXPRESSION:?Set an absolute prepared expression matrix path}"
: "${MANIFEST:?Set an absolute prepared sample manifest path}"
: "${OUTPUT:?Set a unique absolute smoke receipt path}"

cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
# The smoke deliberately uses HiGHS so it neither opens nor depends on a
# Gurobi WLS session while the long metabolic jobs are active.
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME

test -x "${PY}"
test -s "${EXPRESSION}"
test -s "${MANIFEST}"
test -s "${SOURCE_DIR}/collectri_tf_target.tsv.gz"
test -s "${SOURCE_DIR}/omnipath_signed_pkn.tsv.gz"
test ! -e "${OUTPUT}"
mkdir -p "$(dirname "${OUTPUT}")"

srun "${PY}" scripts/run_corneto_regulatory_pilot.py \
    --expression "${EXPRESSION}" \
    --manifest "${MANIFEST}" \
    --collectri "${SOURCE_DIR}/collectri_tf_target.tsv.gz" \
    --pkn "${SOURCE_DIR}/omnipath_signed_pkn.tsv.gz" \
    --study "${STUDY_ACCESSION}" \
    --output "${OUTPUT}" \
    --max-samples "${MAX_SAMPLES:-2}" \
    --max-inputs "${MAX_INPUTS:-3}" \
    --max-outputs "${MAX_OUTPUTS:-6}" \
    --max-depth "${MAX_DEPTH:-3}" \
    --max-edges "${MAX_EDGES:-3000}" \
    --max-flow "${MAX_FLOW:-100}" \
    --lambda-reg "${LAMBDA_REG:-0}" \
    --max-seconds "${MAX_SECONDS:-300}" \
    --solver highs

"${PY}" -c '
import json, pathlib, sys
p = pathlib.Path(sys.argv[1])
r = json.loads(p.read_text())
assert r["study"] == sys.argv[2]
assert r["solver"]["requested"] == "highs"
assert r["solver"]["selected"] == "HIGHS"
assert r["status"] == "completed"
assert 1 <= len(r["samples"]) <= int(sys.argv[3])
assert any(s["status"] in {"optimal", "optimal_inaccurate"} for s in r["samples"])
' "${OUTPUT}" "${STUDY_ACCESSION}" "${MAX_SAMPLES:-2}"
