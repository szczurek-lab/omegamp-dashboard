#!/usr/bin/env bash
# Score the figure 3 flavors FASTAs with APEX, one .tsv per cell.
# Run after gen_figure3_flavors.sh. Needs the BattleAMP environment.
# Scores ALL flavor FASTAs found under flavors/ -- the OmegAMP flavors AND the
# HydrAMP/Prototype baselines -- so panels C/D are compared on one scorer (APEX).
# (Baseline FASTAs must be present; their MIC is loaded by figure3_flavors.py.)
# Usage: bash scripts/reproduce/score_figure3_flavors.sh [--only flavor]
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

flavors_dir="${DATA_DIR}/figure3/flavors"

for fasta in $(find "${flavors_dir}" -name '*.fasta' -not -path '*/_inputs/*' | sort); do
    flavor=$(basename "$(dirname "${fasta}")")
    [ -n "${only}" ] && [ "${only}" != "${flavor}" ] && continue

    stem=$(basename "${fasta}" .fasta)
    pred_dir="$(dirname "${fasta}")/predictions"
    mkdir -p "${pred_dir}"
    prediction="${pred_dir}/${stem}.tsv"
    [ -s "${prediction}" ] && continue

    dataset="figure3-flavors-${flavor}-${stem}"
    cp "${fasta}" "${BATTLEAMP_DIR}/datasets/${dataset}.fasta"

    cd "${BATTLEAMP_DIR}"
    CUDA_VISIBLE_DEVICES=${GPU} snakemake --profile profile/ score \
        --config fasta="datasets/${dataset}.fasta" run_models=apex

    cp "results/inference/apex-min/${dataset}/predictions.tsv" "${prediction}"
done

echo "Done. $(find "${flavors_dir}" -name '*.tsv' | wc -l) prediction files in ${flavors_dir}"
