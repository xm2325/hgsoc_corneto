#!/usr/bin/env bash
# Add the organoid-model comparison to the already validated paired-tissue and
# pseudobulk groups; output is a new v2 receipt and never overwrites v1.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --job-name=ext-reg-compare-v2
#SBATCH --output=/scratch/project_2012997/xiaomei/hgsoc_corneto/logs/ext-reg-compare-v2-%j.out

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=${PYTHON_BIN:-/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python}
EXT=${EXTERNAL_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto_external}
SIG=${ROOT}/data/processed/corneto/external_validation/taylor_regulatory_signature_l0p001.tsv
G277=${EXT}/GSE277107
G189=${EXT}/GSE189955
G208=${EXT}/GSE208216

cd "${ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src"
unset GRB_LICENSE_FILE GRB_LICENSE_SERVER GUROBI_HOME
test -s "${EXT}/regulatory_external_comparison_v1.json"
test ! -e "${G208}/corneto/evidence_v2"
test ! -e "${EXT}/regulatory_external_comparison_v2.json"

srun "${PY}" scripts/build_external_regulatory_evidence.py \
    --solution "${G208}/corneto/joint_l0p001_gurobi_eval_v2.json" \
    --bundle "${G208}/corneto/joint_input_bundle_eval_v2.json" \
    --signature "${SIG}" --manifest "${G208}/corneto_input_v1/corneto_manifest.tsv" \
    --normalization-receipt "${G208}/corneto_input_v1/receipt.json" \
    --study GSE208216 --patient-id-field patient_id --analysis-unit patient_model \
    --group gse208216_pdo=comparison_role:hgsoc_organoid \
    --group gse208216_ft=comparison_role:fallopian_tube_organoid \
    --within-patient-threshold 0.5 --output-dir "${G208}/corneto/evidence_v2"

srun "${PY}" scripts/compare_external_corneto.py \
    --signature "${SIG}" \
    --group "gse277107_ovary=${G277}/corneto/evidence_v2/gse277107_ovary.evidence.tsv,${G277}/corneto/evidence_v2/gse277107_ovary.contract.json" \
    --group "gse277107_omentum=${G277}/corneto/evidence_v2/gse277107_omentum.evidence.tsv,${G277}/corneto/evidence_v2/gse277107_omentum.contract.json" \
    --group "gse189955_hgsoc_epithelial=${G189}/corneto/evidence_v2/gse189955_hgsoc_epithelial.evidence.tsv,${G189}/corneto/evidence_v2/gse189955_hgsoc_epithelial.contract.json" \
    --group "gse189955_fibroblast=${G189}/corneto/evidence_v2/gse189955_fibroblast.evidence.tsv,${G189}/corneto/evidence_v2/gse189955_fibroblast.contract.json" \
    --group "gse189955_normal_ft=${G189}/corneto/evidence_v2/gse189955_normal_ft.evidence.tsv,${G189}/corneto/evidence_v2/gse189955_normal_ft.contract.json" \
    --group "gse208216_pdo=${G208}/corneto/evidence_v2/gse208216_pdo.evidence.tsv,${G208}/corneto/evidence_v2/gse208216_pdo.contract.json" \
    --group "gse208216_ft=${G208}/corneto/evidence_v2/gse208216_ft.evidence.tsv,${G208}/corneto/evidence_v2/gse208216_ft.contract.json" \
    --bootstrap-iterations 2000 --seed 1729 --prevalence-threshold 0.5 \
    --min-patients 2 --output "${EXT}/regulatory_external_comparison_v2.json"
