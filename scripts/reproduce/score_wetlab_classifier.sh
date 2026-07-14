#!/usr/bin/env bash
# Classifier ensemble scoring (broad AMP + species + strains) for wet-lab
# cells. One predictions.csv per cell with the full per-sequence ensemble
# output, no top-K filtering. Needs the OmegAMP environment.
# Usage: bash scripts/reproduce/score_wetlab_classifier.sh [--only motif|analog|denovo]
set -euo pipefail

source config.sh

PREDICT_SCRIPT="${OMEGAMP_DIR}/project/scripts/inference/predict_sequences.py"

only=""
while [ $# -gt 0 ]; do
    case "$1" in
        --only) only="$2"; shift ;;
    esac
    shift
done

wetlab="${DATA_DIR}/wetlab"

for fasta in $(find "${wetlab}" -name 'samples.fasta' | sort); do
    cell="$(dirname "${fasta}")"
    campaign=$(echo "${cell}" | sed -n 's|.*/wetlab/\([^/]*\)/.*|\1|p')
    [ -n "${only}" ] && [ "${only}" != "${campaign}" ] && continue

    predictions="${cell}/predictions.csv"
    [ -s "${predictions}" ] && continue

    cd "${OMEGAMP_DIR}"
    CUDA_VISIBLE_DEVICES=${GPU} python "${PREDICT_SCRIPT}" \
        "${fasta}" \
        --classifier "all" \
        --output_csv "${predictions}"
done

echo "Done. $(find "${wetlab}" -name 'predictions.csv' | wc -l) classifier prediction files in ${wetlab}"
