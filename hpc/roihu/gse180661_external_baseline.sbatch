#!/bin/bash
# Three non-Gurobi stages for GSE180661. Submit each stage explicitly; this
# file does not poll Slurm and does not submit successors by itself.
#
# Stage 1 (small public metadata only):
#   sbatch --export=ALL,STAGE=metadata,MIN_CELLS=20 \
#     hpc/roihu/gse180661_external_baseline.sbatch
# Stage 2 (opt-in 32.6-GB transfer + SHA observation; depend on stage 1):
#   sbatch --dependency=afterok:<metadata_job> --export=ALL,STAGE=matrix_fetch \
#     hpc/roihu/gse180661_external_baseline.sbatch
# Review and commit the generated frozen manifest before Stage 3. This manual
# provenance gate is intentional: the upstream GEO H5 has no published digest.
# Stage 3 (H5 schema/barcode gate + pseudobulk; use the frozen manifest):
#   sbatch --dependency=afterok:<matrix_fetch_job> \
#     --export=ALL,STAGE=pseudobulk,MIN_CELLS=20,FROZEN_MANIFEST=/path/to/frozen.json \
#     hpc/roihu/gse180661_external_baseline.sbatch
#SBATCH --job-name=gse180661-ext
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/gse180661-external-%A.out
#SBATCH --error=logs/gse180661-external-%A.err

set -euo pipefail

: "${PROJECT_ROOT:?Set PROJECT_ROOT to the deployed hgsoc_corneto checkout}"
: "${STAGE:?Set STAGE to metadata, matrix_fetch, or pseudobulk}"

cd "${PROJECT_ROOT}"
mkdir -p logs

PYTHON_BIN="${PYTHON_BIN:-${PROJECT_ROOT}/.venv/bin/python}"
RAW_DIR="${RAW_DIR:-${PROJECT_ROOT}/data/raw/external/GSE180661}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_ROOT}/data/processed/external/GSE180661}"
SOURCE_MANIFEST="${SOURCE_MANIFEST:-${PROJECT_ROOT}/config/external/gse180661_sources.json}"

case "${STAGE}" in
  metadata)
    "${PYTHON_BIN}" scripts/fetch_gse180661_sources.py \
      --manifest "${SOURCE_MANIFEST}" \
      --output-dir "${RAW_DIR}"
    "${PYTHON_BIN}" scripts/prepare_gse180661_external_baseline.py \
      --manifest "${SOURCE_MANIFEST}" \
      --raw-dir "${RAW_DIR}" \
      --output-dir "${OUTPUT_DIR}" \
      --min-cells "${MIN_CELLS:?Set MIN_CELLS before submission}"
    ;;
  matrix_fetch)
    FROZEN_OUTPUT="${FROZEN_OUTPUT:-${OUTPUT_DIR}/gse180661_sources.matrix_sha_observed.json}"
    "${PYTHON_BIN}" scripts/fetch_gse180661_sources.py \
      --manifest "${SOURCE_MANIFEST}" \
      --output-dir "${RAW_DIR}" \
      --include-matrix \
      --write-frozen-manifest "${FROZEN_OUTPUT}"
    echo "Review and commit the observed-SHA manifest before pseudobulk: ${FROZEN_OUTPUT}"
    ;;
  pseudobulk)
    : "${FROZEN_MANIFEST:?Set FROZEN_MANIFEST to the reviewed manifest with H5 SHA-256}"
    "${PYTHON_BIN}" scripts/prepare_gse180661_external_baseline.py \
      --manifest "${FROZEN_MANIFEST}" \
      --raw-dir "${RAW_DIR}" \
      --output-dir "${OUTPUT_DIR}" \
      --min-cells "${MIN_CELLS:?Set MIN_CELLS before submission}" \
      --aggregate \
      --chunk-cells "${CHUNK_CELLS:-4096}"
    ;;
  *)
    echo "Unknown STAGE=${STAGE}; expected metadata, matrix_fetch, or pseudobulk" >&2
    exit 2
    ;;
esac
