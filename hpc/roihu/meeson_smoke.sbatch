#!/usr/bin/env bash
#SBATCH --job-name=hgsoc-meeson-smoke
#SBATCH --account=project_2012997
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4000M
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/meeson-smoke-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
INPUT_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs
VENV_ROOT=${REPO_ROOT}/.venv-cpu
MEESON_COMMIT=1b0e9050720be8fb6769a9a7ab14a4b19bbce8c3
MEESON_CHECKOUT=${INPUT_ROOT}/sources/PhD_2024-${MEESON_COMMIT}
HUMAN_GEM=${INPUT_ROOT}/downloads/Human-GEM-v1.4.1.xml

if ! type module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/zz-csc-env.sh
    set -u
fi
module load python-data/3.12
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1
cd "${REPO_ROOT}"

srun "${VENV_ROOT}/bin/python" -m pytest tests -q
srun "${VENV_ROOT}/bin/python" scripts/audit_meeson_upstream.py \
    --checkout "${MEESON_CHECKOUT}" \
    --human-gem "${HUMAN_GEM}" \
    --expected-commit "${MEESON_COMMIT}" \
    --output-dir data/processed/meeson
srun "${VENV_ROOT}/bin/python" scripts/run_toy_order_benchmark.py \
    --semantics bounds_safe \
    --output data/processed/meeson/toy_order_bounds_safe.json
srun "${VENV_ROOT}/bin/python" scripts/run_toy_order_benchmark.py \
    --semantics published \
    --output data/processed/meeson/toy_order_published.json

printf 'completed_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    > data/processed/meeson/roihu_smoke_completion.tsv
