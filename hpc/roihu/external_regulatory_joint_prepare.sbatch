#!/usr/bin/env bash
# Required exports: STUDY_ACCESSION, EXPRESSION, MANIFEST, OUTPUT_BUNDLE,
# EXPECTED_COUNT, PATIENT_ID_FIELD. Optional INCLUDE_ROLES is comma-separated.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=ext-reg-prepare
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-reg-prepare-%j.out

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
SOURCE_DIR=${SOURCE_DIR:-/scratch/project_2012997/xiaomei/hgsoc_corneto_regulatory}
: "${STUDY_ACCESSION:?}"
: "${EXPRESSION:?}"
: "${MANIFEST:?}"
: "${OUTPUT_BUNDLE:?}"
: "${EXPECTED_COUNT:?}"
: "${PATIENT_ID_FIELD:?}"

cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME
test ! -e "${OUTPUT_BUNDLE}"
mkdir -p "$(dirname "${OUTPUT_BUNDLE}")"

args=(
    --expression "${EXPRESSION}"
    --manifest "${MANIFEST}"
    --collectri "${SOURCE_DIR}/collectri_tf_target.tsv.gz"
    --pkn "${SOURCE_DIR}/omnipath_signed_pkn.tsv.gz"
    --study "${STUDY_ACCESSION}"
    --output "${OUTPUT_BUNDLE}"
    --expected-count "${EXPECTED_COUNT}"
    --patient-id-field "${PATIENT_ID_FIELD}"
    --min-targets 5 --max-inputs 3 --max-outputs 6 --max-depth 3 --max-edges 3000
)
if [[ -n "${INCLUDE_ROLES:-}" ]]; then
    IFS=',' read -r -a roles <<< "${INCLUDE_ROLES}"
    for role in "${roles[@]}"; do
        args+=(--include-role "${role}")
    done
fi
srun "${PY}" scripts/prepare_external_regulatory_multisample.py "${args[@]}"

"${PY}" -c '
import json, sys
r=json.load(open(sys.argv[1]))
assert r["status"] == "completed" and r["response_blind"] is True
assert r["external_validation"] is True
assert r["input_counts"]["selected_conditions"] == int(sys.argv[2])
assert r["input_counts"]["blocked_conditions"] == 0
' "${OUTPUT_BUNDLE}" "${EXPECTED_COUNT}"
