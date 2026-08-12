#!/usr/bin/env bash
# Build patient-level evidence from the two formal external joint fits, then
# compare all five predeclared groups against the frozen Taylor signature.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --job-name=ext-reg-compare
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-reg-compare-%j.out

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
EXT=${EXTERNAL_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}
SIG=${ROOT}/data/processed/corneto/external_validation/taylor_regulatory_signature_l0p001.tsv
G277=${EXT}/GSE277107
G189=${EXT}/GSE189955

cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME
test ! -e "${EXT}/regulatory_external_comparison_v1.json"
test ! -e "${G277}/corneto/evidence_v2"
test ! -e "${G189}/corneto/evidence_v2"

srun "${PY}" scripts/build_external_regulatory_evidence.py \
    --solution "${G277}/corneto/joint_l0p001_gurobi_eval_v2.json" \
    --bundle "${G277}/corneto/joint_input_bundle_eval_v2.json" \
    --signature "${SIG}" \
    --manifest "${G277}/processed/v1/paired_sample_manifest.tsv" \
    --normalization-receipt "${G277}/processed/v1/receipt.json" \
    --study GSE277107 --patient-id-field pair_id --analysis-unit patient_tissue \
    --group gse277107_ovary=normalized_site:ovary \
    --group gse277107_omentum=normalized_site:omentum \
    --within-patient-threshold 0.5 --output-dir "${G277}/corneto/evidence_v2"

srun "${PY}" scripts/build_external_regulatory_evidence.py \
    --solution "${G189}/corneto/joint_l0p001_gurobi_epithelial_fibroblast_eval_v2.json" \
    --bundle "${G189}/corneto/joint_input_bundle_epithelial_fibroblast_eval_v2.json" \
    --signature "${SIG}" \
    --manifest "${G189}/processed/corneto_manifest.tsv" \
    --normalization-receipt "${G189}/processed/patient_celltype_site_cpm_receipt.json" \
    --study GSE189955 --patient-id-field patient_id --analysis-unit patient_pseudobulk \
    --group gse189955_hgsoc_epithelial=comparison_role:hgsoc_epithelial_candidate \
    --group gse189955_fibroblast=comparison_role:hgsoc_site_fibroblast_reference \
    --group gse189955_normal_ft=comparison_role:normal_ft_epithelial_reference \
    --within-patient-threshold 0.5 --output-dir "${G189}/corneto/evidence_v2"

srun "${PY}" scripts/compare_external_corneto.py \
    --signature "${SIG}" \
    --group "gse277107_ovary=${G277}/corneto/evidence_v2/gse277107_ovary.evidence.tsv,${G277}/corneto/evidence_v2/gse277107_ovary.contract.json" \
    --group "gse277107_omentum=${G277}/corneto/evidence_v2/gse277107_omentum.evidence.tsv,${G277}/corneto/evidence_v2/gse277107_omentum.contract.json" \
    --group "gse189955_hgsoc_epithelial=${G189}/corneto/evidence_v2/gse189955_hgsoc_epithelial.evidence.tsv,${G189}/corneto/evidence_v2/gse189955_hgsoc_epithelial.contract.json" \
    --group "gse189955_fibroblast=${G189}/corneto/evidence_v2/gse189955_fibroblast.evidence.tsv,${G189}/corneto/evidence_v2/gse189955_fibroblast.contract.json" \
    --group "gse189955_normal_ft=${G189}/corneto/evidence_v2/gse189955_normal_ft.evidence.tsv,${G189}/corneto/evidence_v2/gse189955_normal_ft.contract.json" \
    --bootstrap-iterations 2000 --seed 1729 --prevalence-threshold 0.5 \
    --min-patients 2 --output "${EXT}/regulatory_external_comparison_v1.json"

"${PY}" -c '
import json, sys
r=json.load(open(sys.argv[1]))
assert r["status"] == "completed"
assert set(r["groups"]) == {
  "gse277107_ovary", "gse277107_omentum", "gse189955_hgsoc_epithelial",
  "gse189955_fibroblast", "gse189955_normal_ft"
}
assert all(g["patient_count"] >= 2 for g in r["groups"].values())
' "${EXT}/regulatory_external_comparison_v1.json"
