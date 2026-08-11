#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=${1:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
INPUT_ROOT=${2:-/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs}
VENV_ROOT=${VENV_ROOT:-${REPO_ROOT}/.venv-cpu}
CORNETO_COMMIT=c2d24f6c914e2d9fd3ebd5a19fc566f9ddc180a8

if ! type module >/dev/null 2>&1; then
    # Non-interactive SSH commands do not always initialize Lmod.
    set +u
    source /etc/profile
    set -u
fi
module load python-data/3.12

python3 -m venv --system-site-packages "${VENV_ROOT}"
"${VENV_ROOT}/bin/python" -m pip install --upgrade pip setuptools wheel
"${VENV_ROOT}/bin/python" -m pip install -e "${REPO_ROOT}[dev,metadata,metabolic]"
"${VENV_ROOT}/bin/python" -m pip install \
    "corneto @ git+https://github.com/saezlab/corneto@${CORNETO_COMMIT}" \
    highspy

"${VENV_ROOT}/bin/python" "${REPO_ROOT}/scripts/fetch_metabolic_sources.py" \
    --destination "${INPUT_ROOT}"

mkdir -p "${REPO_ROOT}/logs" "${REPO_ROOT}/data/processed/meeson"
{
    printf 'repo_commit\t'
    git -C "${REPO_ROOT}" rev-parse HEAD
    printf 'python\t'
    "${VENV_ROOT}/bin/python" --version
    "${VENV_ROOT}/bin/python" -c \
        'from importlib.metadata import version; import highspy; print("cobra\t" + version("cobra")); print("corneto\t" + version("corneto")); print("highspy\tinstalled")'
    printf 'human_gem_sha256\t'
    sha256sum "${INPUT_ROOT}/downloads/Human-GEM-v1.4.1.xml" | cut -d' ' -f1
} > "${REPO_ROOT}/data/processed/meeson/roihu_cpu_environment.tsv"

"${VENV_ROOT}/bin/python" -m pytest "${REPO_ROOT}/tests" -q
