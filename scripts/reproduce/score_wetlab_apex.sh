#!/usr/bin/env bash
# APEX scoring for wet-lab cells. One apex.tsv per cell, with per-strain
# columns plus the min. Needs the BattleAMP environment.
# Usage: bash scripts/reproduce/score_wetlab_apex.sh [--only motif|analog|denovo]
set -euo pipefail

source config.sh
source "${BATTLEAMP_ACTIVATE}"

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

    prediction="${cell}/apex.tsv"
    [ -s "${prediction}" ] && continue

    # BattleAMP outputs to results/inference/apex/<dataset>/; give each cell a
    # unique dataset name when staging.
    rel=$(realpath --relative-to="${wetlab}" "${cell}")
    dataset="wetlab-$(echo "${rel}" | tr '/' '-')"
    cp "${fasta}" "${BATTLEAMP_DIR}/datasets/${dataset}.fasta"

    cd "${BATTLEAMP_DIR}"
    CUDA_VISIBLE_DEVICES=${GPU} snakemake --profile profile/ score \
        --config fasta="datasets/${dataset}.fasta" run_models=apex

    cp "results/inference/apex/${dataset}/predictions.tsv" "${prediction}"
done

echo "Done. $(find "${wetlab}" -name 'apex.tsv' | wc -l) APEX prediction files in ${wetlab}"
