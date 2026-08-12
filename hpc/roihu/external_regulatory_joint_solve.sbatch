#!/usr/bin/env bash
# Required exports: BUNDLE, OUTPUT. Formal external fit at the frozen Taylor
# nominal lambda, independently normalized within the external dataset.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=ext-reg-joint
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-reg-joint-%j.out

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
: "${BUNDLE:?}"
: "${OUTPUT:?}"
cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
export GRB_LICENSE_FILE=${GRB_LICENSE_FILE:-${ROOT}/gurobi_xm.lic}
test ! -e "${OUTPUT}"
mkdir -p "$(dirname "${OUTPUT}")"

srun "${PY}" scripts/run_corneto_regulatory_multisample.py \
    --bundle "${BUNDLE}" --output "${OUTPUT}" --mode pooled \
    --lambda-reg 0.001 --lambda-scaling mean_fit --solver gurobi \
    --max-flow 100 --max-seconds 82800

"${PY}" -c '
import json, sys
r=json.load(open(sys.argv[1]))
assert r["status"] == "completed"
assert r["response_blind"] is True and r["method"]["single_joint_problem"] is True
assert r["method"]["lambda_nominal"] == 0.001
assert r["method"]["lambda_scaling"] == "mean_fit"
assert r["scope_counts"]["preprocessing_blocked"] == 0
assert r["solver"]["selected"] == "GUROBI" and r["solver"]["has_incumbent"] is True
assert r["solver"]["status"] in {"optimal", "optimal_inaccurate"}
assert all(c["status"] in {"optimal", "optimal_inaccurate"} for c in r["conditions"])
' "${OUTPUT}"
