#!/usr/bin/env bash
# Final retry after normalizing the current NMF ``state`` column.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=reg-nmf-all-v4

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src"
OUTROOT=data/processed/corneto/regulatory_nmf_integration_all_cohorts/l0p001_v4
test ! -e "$OUTROOT"
mkdir -p "$OUTROOT"
STUDIES=(E-MTAB-7223 E-MTAB-10801 E-MTAB-11000 E-MTAB-14568)
for study in "${STUDIES[@]}"; do
  receipt=data/processed/corneto/regulatory_multisample/v1/grid/$study/l0p001.json
  assignments=data/processed/rna/$study/nmf_primary/rank_3/assignments.tsv
  test -s "$receipt"; test -s "$assignments"
  srun --exclusive "$PY" scripts/integrate_regulatory_nmf_states.py \
    --regulatory "$receipt" --assignments "$assignments" --study "$study" \
    --output "$OUTROOT/${study}.json"
done
exec "$PY" -c 'import json; from pathlib import Path
root=Path("data/processed/corneto/regulatory_nmf_integration_all_cohorts/l0p001_v4")
rows=[]
for p in sorted(root.glob("E-MTAB-*.json")):
    d=json.loads(p.read_text()); rows.append({"study":d["study_accession"],"matched_samples":d["coverage"]["matched_samples"],"states":len(d["state_summaries"]),"edge_state_tests":len(d["edge_state_tests"]),"path":str(p)})
(root/"summary.json").write_text(json.dumps({"status":"completed","schema_version":"regulatory_nmf_all_cohorts.v1","response_blind":True,"lambda_nominal":0.001,"cohorts":rows,"claim_limit":"descriptive state-edge association only; no phenotype or causal interpretation"},indent=2,sort_keys=True)+"\n")
print(json.dumps({"status":"completed","cohort_count":len(rows),"matched_samples":sum(r["matched_samples"] for r in rows)}))'
