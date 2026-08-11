#!/usr/bin/env bash
#SBATCH --job-name=tpi1-fva-preflight
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/tpi1-fva-preflight-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
INPUT_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs
PYTHON_BIN=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python

cd "${REPO_ROOT}"
srun "${PYTHON_BIN}" scripts/prepare_tpi1_fva_preflight.py \
    --sbml "${INPUT_ROOT}/downloads/Human-GEM-v1.4.1.xml" \
    --symbol TPI1 \
    --biomass-reaction biomass_human \
    --receipt E-MTAB-7223=data/processed/corneto/E-MTAB-7223/full_direct_b25.json \
    --receipt E-MTAB-10801=data/processed/corneto/E-MTAB-10801/full_direct_b25.json \
    --receipt E-MTAB-11000=data/processed/corneto/E-MTAB-11000/full_direct_b25.json \
    --receipt E-MTAB-14568=data/processed/corneto/E-MTAB-14568/full_direct_b25.json \
    --output data/processed/corneto/tpi1_fva_preflight.json
