#!/usr/bin/env bash
#SBATCH --job-name=tpi1-fba-fva
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/tpi1-fba-fva-%j.out

set -euo pipefail
ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
INPUT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
cd "${ROOT}"
exec srun "${PY}" scripts/run_tpi1_fba_fva.py \
  --preflight data/processed/corneto/tpi1_fva_preflight.json \
  --sbml "${INPUT}/downloads/Human-GEM-v1.4.1.xml" \
  --output data/processed/meeson/tpi1_human_gem_fba_fva.json
