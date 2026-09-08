#!/usr/bin/env bash
# All modes are licence-free and isolated from the frozen b25 recovery chain.
#SBATCH --account=project_2012997
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --job-name=post-audit-test

set -euo pipefail
ROOT=${REPO_ROOT:-/scratch/project_2012997/xiaomei/hgsoc_corneto}
PY=/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python
RNA=/scratch/project_2012997/xiaomei/hgsoc_corneto_rna/aggregated
MODEL=/scratch/project_2012997/xiaomei/hgsoc_corneto_inputs/downloads/Human-GEM-v1.4.1.xml
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONPATH="$ROOT/src:$ROOT/scripts"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
MODE=${VALIDATION_MODE:?Set VALIDATION_MODE}
VALIDATION_TAG=${VALIDATION_TAG:-}
[[ "$VALIDATION_TAG" =~ ^[a-zA-Z0-9_-]*$ ]]
OUT="data/processed/post_audit_validation_20260908/${MODE}${VALIDATION_TAG}"
test ! -e "$OUT"
case "$MODE" in
  nmf_smoke)
    exec srun "$PY" scripts/run_patient_grouped_nmf_validation.py --root "$ROOT" --rna-root "$RNA" \
      --output-dir "$OUT" --top-genes 500 --starts 2 --draws 2 ;;
  nmf_validation)
    "$PY" -c 'import json; from pathlib import Path; p=Path("data/processed/post_audit_validation_20260908/nmf_smoke/receipt.json"); d=json.loads(p.read_text()); assert d["status"]=="completed" and len(d["folds"])==8'
    exec srun "$PY" scripts/run_patient_grouped_nmf_validation.py --root "$ROOT" --rna-root "$RNA" \
      --output-dir "$OUT" --top-genes 6000 --starts 20 --draws 30 ;;
  metabolic_smoke)
    exec srun "$PY" scripts/run_metabolic_information_pilot.py --root "$ROOT" --rna-root "$RNA" \
      --human-gem "$MODEL" --output-dir "$OUT" --per-study 1 --shuffles 2 ;;
  metabolic_pilot)
    "$PY" -c 'import json,sys; from pathlib import Path; p=Path(sys.argv[1]); d=json.loads(p.read_text()); assert d["status"]=="completed" and len(d["sample_information"])==4' "data/processed/post_audit_validation_20260908/metabolic_smoke${VALIDATION_TAG}/receipt.json"
    exec srun "$PY" scripts/run_metabolic_information_pilot.py --root "$ROOT" --rna-root "$RNA" \
      --human-gem "$MODEL" --output-dir "$OUT" --per-study 2 --shuffles 10 ;;
  biology)
    GMT="$ROOT/data/inputs/post_audit_validation_20260908/h.all.v2025.1.Hs.symbols.gmt"
    test -s "$GMT"
    "$PY" -c 'import hashlib,sys; assert hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest()=="f22066af72e215ccb7b89d88e492c07e1eef17534c2ca7b0f9902cfecbbdd8e9"' "$GMT"
    exec srun "$PY" scripts/interpret_existing_biology.py --root "$ROOT" --rna-root "$RNA" \
      --hallmark-gmt "$GMT" --output-dir "$OUT" ;;
  *) echo "Unknown validation mode: $MODE" >&2; exit 2 ;;
esac
