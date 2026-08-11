#!/usr/bin/env bash
#SBATCH --job-name=hgsoc-corneto-14568
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/corneto-14568-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
RNA_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna
INPUT_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs
VENV_ROOT=${HGSOC_CORNETO_ENV:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env}
SOLVER=${SOLVER:-auto}
MOSEK_LICENSE=${MOSEKLM_LICENSE_FILE:-${REPO_ROOT}/mosek_xm.lic}
GUROBI_LICENSE=${GRB_LICENSE_FILE:-${REPO_ROOT}/gurobi_xm.lic}
PILOT_OUTPUT=${PILOT_OUTPUT:-data/processed/corneto/E-MTAB-14568/pilot_receipt.json}

if ! type module >/dev/null 2>&1; then
    set +u
    export CSC_ENV_INIT_NON_INTERACTIVE=yes
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
if [[ "${SOLVER}" == "mosek" && ! -r "${MOSEK_LICENSE}" ]]; then
    echo "SOLVER=mosek but MOSEK license is not readable: ${MOSEK_LICENSE}" >&2
    exit 2
fi
if [[ "${SOLVER}" == "gurobi" && ! -r "${GUROBI_LICENSE}" ]]; then
    echo "SOLVER=gurobi but Gurobi license is not readable: ${GUROBI_LICENSE}" >&2
    exit 2
fi
if [[ -r "${MOSEK_LICENSE}" ]]; then
    export MOSEKLM_LICENSE_FILE="${MOSEK_LICENSE}"
fi
if [[ -r "${GUROBI_LICENSE}" ]]; then
    export GRB_LICENSE_FILE="${GUROBI_LICENSE}"
fi
export HGSOC_CORNETO_REPO_COMMIT=$(git -C "${REPO_ROOT}" rev-parse HEAD)
cd "${REPO_ROOT}"

srun "${VENV_ROOT}/bin/python" scripts/run_corneto_14568_pilot.py \
    --expression "${RNA_ROOT}/aggregated/E-MTAB-14568/gene_tpm.tsv.gz" \
    --manifest data/processed/metadata/ocm_master_manifest.tsv \
    --human-gem "${INPUT_ROOT}/downloads/Human-GEM-v1.4.1.xml" \
    --output "${PILOT_OUTPUT}" \
    --primary-only \
    --max-samples "${MAX_SAMPLES:-3}" \
    --max-candidates "${MAX_CANDIDATES:-250}" \
    --expression-transform "${EXPRESSION_TRANSFORM:-log1p_tpm}" \
    --solver "${SOLVER}"
