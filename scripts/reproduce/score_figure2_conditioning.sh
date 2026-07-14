#!/usr/bin/env bash
# APEX MIC predictions for every figure 2 conditioning FASTA.
# Reads each samples.fasta under data/figure2/conditioning/ and writes
# predictions.tsv next to it.
# Usage: bash scripts/reproduce/score_figure2_conditioning.sh [--mini]
set -euo pipefail

source config.sh
source "${BATTLEAMP_ACTIVATE}"

samples_dirs=$(find "${DATA_DIR}/figure2/conditioning" -mindepth 2 -maxdepth 2 -type d | sort)
[ "${1:-}" = "--mini" ] && samples_dirs=$(echo "${samples_dirs}" | head -3)

for dir in ${samples_dirs}; do
    [ -s "${dir}/samples.fasta" ] || continue
    [ -s "${dir}/predictions.tsv" ] && continue

    # BattleAMP keys outputs by FASTA filename stem, so stage each cell under a
    # unique stem to avoid collisions across the 170 "samples.fasta" files.
    rel=${dir#${DATA_DIR}/figure2/conditioning/}
    dataset="figure2-${rel//\//_}"
    staged="${BATTLEAMP_DIR}/datasets/${dataset}.fasta"
    cp "${dir}/samples.fasta" "${staged}"

    cd "${BATTLEAMP_DIR}"
    CUDA_VISIBLE_DEVICES=${GPU} snakemake --profile profile/ score \
        --config fasta="${staged}" run_models=apex

    cp "results/inference/apex-min/${dataset}/predictions.tsv" "${dir}/predictions.tsv"
done

echo "Done. $(find "${DATA_DIR}/figure2/conditioning" -name "predictions.tsv" | wc -l) APEX predictions."