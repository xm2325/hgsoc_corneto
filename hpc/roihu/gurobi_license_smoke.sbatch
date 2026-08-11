#!/usr/bin/env bash
#SBATCH --job-name=hgsoc-gurobi-smoke
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/gurobi-license-smoke-%j.out

set -euo pipefail

REPO_ROOT=/scratch/project_2012997/xiaomei/hgsoc_corneto
VENV_ROOT=${HGSOC_CORNETO_ENV:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env}
LICENSE=${GRB_LICENSE_FILE:-${REPO_ROOT}/gurobi_xm.lic}
if [[ ! -r "${LICENSE}" ]]; then
    echo "Gurobi license is not readable" >&2
    exit 2
fi
export GRB_LICENSE_FILE="${LICENSE}"
"${VENV_ROOT}/bin/python" -c 'import gurobipy as gp; env=gp.Env(empty=True); env.setParam("OutputFlag",0); env.start(); model=gp.Model(env=env); x=model.addVar(lb=0,ub=1); model.setObjective(x,gp.GRB.MAXIMIZE); model.optimize(); print("gurobi_compute_smoke="+("ok" if model.Status==gp.GRB.OPTIMAL else "status_"+str(model.Status))); model.dispose(); env.dispose()'
