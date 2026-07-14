#!/usr/bin/env bash
# Score the figure 3 sweep FASTAs with APEX, one .tsv per tau/sigma cell.
# Run after gen_figure3_sweep.sh. Needs the BattleAMP environment.
# Usage: bash scripts/reproduce/score_figure3_sweep.sh [--mini] [--only positive|negative]
set -euo pipefail

source config.sh
source "${BATTLEAMP_ACTIVATE}"

# Which prototype sets to score, and how many cells from each.
sets="positive negative"
max_cells=0   # 0 means all

while [ $# -gt 0 ]; do
    case "$1" in
        --mini) max_cells=3 ;;
        --only) sets="$2"; shift ;;
    esac
    shift
done

sweep="${DATA_DIR}/figure3/sweep"

for set_name in ${sets}; do
    seq_dir="${sweep}/${set_name}/sequences"
    pred_dir="${sweep}/${set_name}/predictions"
    mkdir -p "${pred_dir}"

    fastas=$(find "${seq_dir}" -name 'tau*.fasta' | sort)
    [ "${max_cells}" -gt 0 ] && fastas=$(echo "${fastas}" | head -"${max_cells}")

    for fasta in ${fastas}; do
        cell=$(basename "${fasta}" .fasta)
        prediction="${pred_dir}/${cell}.tsv"
        [ -s "${prediction}" ] && continue

        # BattleAMP names its output after the input FASTA stem, so give each
        # cell a unique stem when staging it.
        dataset="figure3-sweep-${set_name}-${cell}"
        cp "${fasta}" "${BATTLEAMP_DIR}/datasets/${dataset}.fasta"

        cd "${BATTLEAMP_DIR}"
        CUDA_VISIBLE_DEVICES=${GPU} snakemake --profile profile/ score \
            --config fasta="datasets/${dataset}.fasta" run_models=apex

        cp "results/inference/apex-min/${dataset}/predictions.tsv" "${prediction}"
    done
done

echo "Done. $(find "${sweep}" -name '*.tsv' | wc -l) prediction files in ${sweep}"