#!/usr/bin/env bash
#SBATCH --job-name=compare-pooled-cohort-metabolic
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/compare-pooled-cohort-metabolic-%j.out

set -euo pipefail
ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
cd "${ROOT}"
srun "${PY}" scripts/compare_pooled_cohort_metabolic.py \
    --pooled data/processed/corneto/pooled_primary_60/joint_full_60_b25.json \
    --cohort E-MTAB-7223=data/processed/corneto/E-MTAB-7223/full_direct_b25.json \
    --cohort E-MTAB-10801=data/processed/corneto/E-MTAB-10801/full_direct_b25.json \
    --cohort E-MTAB-11000=data/processed/corneto/E-MTAB-11000/full_direct_b25.json \
    --cohort E-MTAB-14568=data/processed/corneto/E-MTAB-14568/full_direct_b25.json \
    --output data/processed/corneto/pooled_vs_cohort_metabolic_b25.json
