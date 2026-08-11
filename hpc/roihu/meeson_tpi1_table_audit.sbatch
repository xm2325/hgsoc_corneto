#!/usr/bin/env bash
#SBATCH --job-name=meeson-tpi1-audit
#SBATCH --account=project_2012997
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/meeson-tpi1-audit-%j.out

set -euo pipefail

ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
INPUT=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs
PY=${ROOT}/.venv-cpu/bin/python
MEESON=PhD_2024-1b0e9050720be8fb6769a9a7ab14a4b19bbce8c3

cd "${ROOT}"
exec srun "${PY}" scripts/audit_meeson_tpi1_validation.py \
  --input "${INPUT}/sources/${MEESON}/siTPI1_validation/siTPI1_experimental_simulation_comparison.csv" \
  --output data/processed/meeson/tpi1_published_table_audit.json
