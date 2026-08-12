#!/usr/bin/env bash
#SBATCH --job-name=metabolic-availability
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/metabolic-availability-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
PYTHON_BIN=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
OUTPUT_PATH="${OUTPUT_PATH:-data/processed/corneto/audits/metabolic_baseline_availability_${SLURM_JOB_ID}.json}"

cd "${REPO_ROOT}"
set +e
srun "${PYTHON_BIN}" scripts/audit_metabolic_baseline_availability.py \
    --receipt E-MTAB-7223=data/processed/corneto/E-MTAB-7223/full_direct_b25.json \
    --receipt E-MTAB-10801=data/processed/corneto/E-MTAB-10801/full_direct_b25.json \
    --receipt E-MTAB-11000=data/processed/corneto/E-MTAB-11000/full_direct_b25.json \
    --receipt E-MTAB-14568=data/processed/corneto/E-MTAB-14568/full_direct_b25.json \
    --output "${OUTPUT_PATH}"
audit_status=$?
set -e

# Exit 0 for an expected incomplete snapshot; retain nonzero for malformed
# existing receipts.  The JSON always records the actual audit status.
if [[ ${audit_status} -eq 2 ]]; then
    exit 0
fi
exit "${audit_status}"
