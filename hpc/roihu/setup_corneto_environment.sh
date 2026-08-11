#!/usr/bin/env bash
# Build a dedicated, reproducible CORNETO environment separate from the RNA venv.

set -euo pipefail

REPO_ROOT=${1:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
ENV_ROOT=${HGSOC_CORNETO_ENV:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env}
CORNETO_COMMIT=c2d24f6c914e2d9fd3ebd5a19fc566f9ddc180a8

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12

if [[ ! -x "${ENV_ROOT}/bin/python" ]]; then
    python3 -m venv --system-site-packages "${ENV_ROOT}"
fi

PYTHON=${ENV_ROOT}/bin/python
"${PYTHON}" -m pip install --upgrade pip setuptools wheel
"${PYTHON}" -m pip install -e "${REPO_ROOT}[analysis,dev,metadata,metabolic]"
"${PYTHON}" -m pip install \
    "corneto @ git+https://github.com/saezlab/corneto@${CORNETO_COMMIT}" \
    highspy \
    'numpy>=2,<2.5'

# Gurobi is optional because its package and license are site/user state.  If
# the caller explicitly requests it, install only the client package and let
# the solver smoke test validate the license without printing its contents.
if [[ "${INSTALL_GUROBI:-0}" == "1" ]]; then
    "${PYTHON}" -m pip install gurobipy
fi

mkdir -p "${REPO_ROOT}/logs" "${REPO_ROOT}/data/processed/corneto"
{
    printf 'repo_commit\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
    printf 'environment_root\t%s\n' "${ENV_ROOT}"
    printf 'python\t%s\n' "$("${PYTHON}" --version 2>&1)"
    printf 'numpy\t%s\n' "$("${PYTHON}" -c 'from importlib.metadata import version; print(version("numpy"))')"
    printf 'cobra\t%s\n' "$("${PYTHON}" -c 'from importlib.metadata import version; print(version("cobra"))')"
    printf 'corneto\t%s\n' "$("${PYTHON}" -c 'from importlib.metadata import version; print(version("corneto"))')"
    printf 'highspy\t%s\n' "$("${PYTHON}" -c 'from importlib.metadata import version; print(version("highspy"))')"
    if "${PYTHON}" -c 'import gurobipy' >/dev/null 2>&1; then
        printf 'gurobipy\tinstalled\n'
    else
        printf 'gurobipy\tnot_installed\n'
    fi
    if [[ -n "${GRB_LICENSE_FILE:-}" ]]; then
        printf 'grb_license_file\tconfigured\n'
    else
        printf 'grb_license_file\tnot_configured\n'
    fi
} > "${REPO_ROOT}/data/processed/corneto/environment_receipt.tsv"

"${PYTHON}" -m pytest "${REPO_ROOT}/tests" -q
