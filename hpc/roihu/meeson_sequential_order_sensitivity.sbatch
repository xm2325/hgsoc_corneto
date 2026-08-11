#!/usr/bin/env bash
#SBATCH --job-name=meeson-order
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/meeson-order-%x-%j.out

set -euo pipefail

STUDY=${1:?supply study accession}
ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
INPUT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs
RNA=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated
PY=${ROOT}/hgsoc_corneto_env/bin/python

if ! type module >/dev/null 2>&1; then
  set +u; source /etc/profile.d/zz-csc-env.sh; set -u
fi
module load python-data/3.12
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1
cd "${ROOT}"
exec srun "${PY}" scripts/run_meeson_sequential_order_sensitivity.py \
  --study "${STUDY}" \
  --manifest data/processed/metadata/ocm_master_manifest.tsv \
  --expression "${RNA}/${STUDY}/gene_tpm.tsv.gz" \
  --human-gem "${INPUT}/downloads/Human-GEM-v1.4.1.xml" \
  --receipt "data/processed/corneto/${STUDY}/full_direct_b25.json" \
  --output "data/processed/meeson/order_sensitivity/${STUDY}_bounds_safe.json" \
  --permutations 25 --semantics bounds_safe
