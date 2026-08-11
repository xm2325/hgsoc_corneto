#!/usr/bin/env bash
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=reg-lambda-compare

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
cd "$ROOT"
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
BASE=data/processed/corneto/E-MTAB-14568/regulatory_full_gurobi_27_l0.json
test -s "$BASE"
for label in l0p01 l0p05 l0p1 l0p25 l0p5 l1 l2; do
  test -s "data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_${label}.json"
done
exec "$PY" scripts/compare_regulatory_lambda_receipts.py \
  --pilot 0="$BASE" \
  --pilot 0.01=data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_l0p01.json \
  --pilot 0.05=data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_l0p05.json \
  --pilot 0.1=data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_l0p1.json \
  --pilot 0.25=data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_l0p25.json \
  --pilot 0.5=data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_l0p5.json \
  --pilot 1=data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_l1.json \
  --pilot 2=data/processed/corneto/E-MTAB-14568/regulatory_robustness_full_27_l2.json \
  --output data/processed/corneto/E-MTAB-14568/regulatory_lambda_robustness_comparison.json
