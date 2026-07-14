#!/usr/bin/env bash
# Analog tau/sigma sweep for figure 3: generation only.
# Generates 10 analogs per prototype across the tau x sigma grid.
# Score separately with score_figure3_sweep.sh. Run under the OmegAMP env.
# Usage: bash scripts/reproduce/gen_figure3_sweep.sh [--mini] [--only positive|negative]
set -euo pipefail

source config.sh

ANALOGS_PER_PROTOTYPE=10
TAU_VALUES=(0 0.05 0.1 0.15 0.2 0.3 0.5 0.7)
SIGMA_VALUES=(0.0 0.25 0.5 0.75 1.0)
ONLY=both

while [ $# -gt 0 ]; do
    case "$1" in
        --mini)
            ANALOGS_PER_PROTOTYPE=2
            TAU_VALUES=(0 0.15 0.5)
            SIGMA_VALUES=(0.0 1.0)
            ;;
        --only) ONLY="$2"; shift ;;
    esac
    shift
done

OUT="${DATA_DIR}/figure3/sweep"

generate_set() {
    local proto_fasta="$1"
    local set_name="$2"
    local seq_dir="${OUT}/${set_name}/sequences"
    mkdir -p "${seq_dir}"

    # Each prototype repeated, so we get ANALOGS_PER_PROTOTYPE analogs of each.
    local repeated="${seq_dir}/prototypes-repeated.fasta"
    repeat_fasta "${proto_fasta}" "${ANALOGS_PER_PROTOTYPE}" "${repeated}"
    local n=$(grep -c '^>' "${repeated}")

    cd "${OMEGAMP_DIR}"
    for tau in "${TAU_VALUES[@]}"; do
        for sigma in "${SIGMA_VALUES[@]}"; do
            local fasta="${seq_dir}/tau${tau}_sigma${sigma}.fasta"
            [ -s "${fasta}" ] && continue
            CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog prototype-derived \
                --checkpoint_path "${CHECKPOINT}" \
                --analog_sequences "${repeated}" \
                --prototype_sequences "${repeated}" \
                --tau "${tau}" --sigma "${sigma}" \
                --num_samples "${n}" --batch_size 256 \
                --output_fasta "${fasta}" \
                --conditioning_output_path "${seq_dir}/tau${tau}_sigma${sigma}_conditioning.pt"
        done
    done
}

if [ "${ONLY}" != "negative" ]; then
    generate_set "${POSITIVE_PROTOTYPES}" positive
fi
if [ "${ONLY}" != "positive" ]; then
    generate_set "${NEGATIVE_PROTOTYPES}" negative
fi

echo "Done. $(find "${OUT}" -name 'tau*.fasta' | wc -l) FASTA files in ${OUT}"