#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=gse277107-intake
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/gse277107-intake-%j.out

set -euo pipefail
umask 027

ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
EXTERNAL_ROOT=${EXTERNAL_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}
STUDY_ROOT=${EXTERNAL_ROOT}/GSE277107
RAW_DIR=${STUDY_ROOT}/raw
PROCESSED_DIR=${STUDY_ROOT}/processed/v1
FAILURE_DIR=${STUDY_ROOT}/failures
STAGING_ROOT=${STUDY_ROOT}/staging
STAGE_DIR=${STAGING_ROOT}/${SLURM_JOB_ID:-manual}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
SOURCE_MANIFEST=${ROOT}/config/external/gse277107_sources.json
GATE_RECEIPT=${PROCESSED_DIR}/receipt_gate.json

mkdir -p "${RAW_DIR}" "${FAILURE_DIR}" "${PROCESSED_DIR%/*}" "${STAGING_ROOT}"
failure_tmp=$(mktemp "${FAILURE_DIR}/.gse277107-${SLURM_JOB_ID:-manual}.XXXXXX")
preserve_failure() {
    exit_code=$?
    if [[ ${exit_code} -ne 0 ]]; then
        {
            printf 'schema_version\tgse277107_hpc_failure.v1\n'
            printf 'status\tfailed\n'
            printf 'job_id\t%s\n' "${SLURM_JOB_ID:-manual}"
            printf 'exit_code\t%s\n' "${exit_code}"
            printf 'failed_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
            printf 'repo_root\t%s\n' "${ROOT}"
            printf 'raw_dir\t%s\n' "${RAW_DIR}"
            printf 'processed_dir\t%s\n' "${PROCESSED_DIR}"
            printf 'stage_dir\t%s\n' "${STAGE_DIR}"
            printf 'slurm_log\t%s\n' "${ROOT}/logs/gse277107-intake-${SLURM_JOB_ID:-manual}.out"
        } > "${failure_tmp}"
        mv "${failure_tmp}" "${FAILURE_DIR}/gse277107-${SLURM_JOB_ID:-manual}.failure.tsv"
    else
        rm -f "${failure_tmp}"
    fi
}
trap preserve_failure EXIT

cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
unset GRB_LICENSE_FILE GUROBI_HOME || true
test -x "${PY}"
test -f "${SOURCE_MANIFEST}"

# A completed fixed output is immutable on rerun: audit it, do not recompute.
if [[ -f "${PROCESSED_DIR}/receipt.json" ]]; then
    srun "${PY}" scripts/audit_gse277107_receipt.py \
        --processed-dir "${PROCESSED_DIR}" \
        --output "${GATE_RECEIPT}"
    exit 0
fi
if [[ -e "${PROCESSED_DIR}" ]]; then
    printf 'Refusing to overwrite incomplete fixed output: %s\n' "${PROCESSED_DIR}" >&2
    exit 67
fi
test ! -e "${STAGE_DIR}"
mkdir "${STAGE_DIR}"

srun "${PY}" scripts/fetch_gse277107.py \
    --source-manifest "${SOURCE_MANIFEST}" \
    --output-dir "${RAW_DIR}"
srun "${PY}" scripts/prepare_gse277107.py \
    --source-manifest "${SOURCE_MANIFEST}" \
    --source-dir "${RAW_DIR}" \
    --output-dir "${STAGE_DIR}"
srun "${PY}" scripts/audit_gse277107_receipt.py \
    --processed-dir "${STAGE_DIR}" \
    --output "${STAGE_DIR}/receipt_gate.json"

# The staging and fixed directories share STUDY_ROOT, so this publishes the
# completely audited bundle with one directory rename.
mv "${STAGE_DIR}" "${PROCESSED_DIR}"
