#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=${1:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
TOOL_ROOT=${2:-/scratch/project_2012997/xiaomei/hgsoc_corneto_tools}
VENV_ROOT=${VENV_ROOT:-${REPO_ROOT}/.venv-cpu}
SALMON_VERSION=1.10.0
SALMON_URL=https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
SALMON_SHA256=b876d041ef3bfbe44422b052b99ce387ff4e521c76002355c7b27882cf19c01b
SALMON_BYTES=93320991
ARCHIVE=${TOOL_ROOT}/downloads/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
SALMON_ROOT=${TOOL_ROOT}/salmon-${SALMON_VERSION}

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12

mkdir -p "${TOOL_ROOT}/downloads" "${REPO_ROOT}/data/processed/rna"

if [[ -f "${ARCHIVE}" ]]; then
    [[ "$(stat -c %s "${ARCHIVE}")" == "${SALMON_BYTES}" ]]
    [[ "$(sha256sum "${ARCHIVE}" | cut -d' ' -f1)" == "${SALMON_SHA256}" ]]
else
    curl --location --fail --retry 5 --retry-delay 10 --continue-at - \
        --output "${ARCHIVE}.partial" "${SALMON_URL}"
    [[ "$(stat -c %s "${ARCHIVE}.partial")" == "${SALMON_BYTES}" ]]
    [[ "$(sha256sum "${ARCHIVE}.partial" | cut -d' ' -f1)" == "${SALMON_SHA256}" ]]
    mv "${ARCHIVE}.partial" "${ARCHIVE}"
fi

if [[ ! -x "${SALMON_ROOT}/bin/salmon" ]]; then
    STAGING=$(mktemp -d "${TOOL_ROOT}/.salmon-${SALMON_VERSION}.XXXXXX")
    trap 'rm -rf "${STAGING}"' EXIT
    tar -xzf "${ARCHIVE}" -C "${STAGING}"
    "${STAGING}/salmon-latest_linux_x86_64/bin/salmon" --version
    mv "${STAGING}/salmon-latest_linux_x86_64" "${SALMON_ROOT}"
    rmdir "${STAGING}"
    trap - EXIT
fi

[[ "$("${SALMON_ROOT}/bin/salmon" --version)" == "salmon ${SALMON_VERSION}" ]]
"${VENV_ROOT}/bin/python" -m pip install -e \
    "${REPO_ROOT}[analysis,dev,metadata,metabolic]"
"${VENV_ROOT}/bin/python" -m pip check

{
    printf 'repo_commit\t'
    git -C "${REPO_ROOT}" rev-parse HEAD
    printf 'python\t'
    "${VENV_ROOT}/bin/python" --version
    "${VENV_ROOT}/bin/python" -c \
        'from importlib.metadata import version; print("numpy\t" + version("numpy")); print("pandas\t" + version("pandas")); print("scikit_learn\t" + version("scikit-learn"))'
    printf 'salmon\t%s\n' "$("${SALMON_ROOT}/bin/salmon" --version)"
    printf 'salmon_archive_bytes\t%s\n' "${SALMON_BYTES}"
    printf 'salmon_archive_sha256\t%s\n' "${SALMON_SHA256}"
} > "${REPO_ROOT}/data/processed/rna/roihu_rna_environment.tsv"
