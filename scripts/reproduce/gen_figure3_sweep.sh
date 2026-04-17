#!/usr/bin/env bash
# =============================================================================
# Generate the analog parameter sweep dataset for Figure 3.
#
# Runs all 8 × 5 = 40 (τ, σ) combinations for both positive and negative
# prototype sets, then runs APEX MIC predictions on each.
#
# Writes:
#   data/figure3/sweep/{positive,negative}/sequences/{pt}_tau{T}_sigma{S}.fasta
#   data/figure3/sweep/{positive,negative}/predictions/{pt}_tau{T}_sigma{S}-min.tsv
#
# Prototype files (50 sequences each by default):
#   OMEGAMP_DIR/data/activity-data/curated-AMPs_samples.fasta
#   OMEGAMP_DIR/data/activity-data/curated-Non-AMPs_samples.fasta
#
# To use the full curated datasets instead of the samples, set:
#   POSITIVE_PROTOTYPES=/path/to/curated-AMPs.fasta
#   NEGATIVE_PROTOTYPES=/path/to/curated-Non-AMPs.fasta
#
# Runtime: several hours on a single GPU for the full τ range.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/config.sh"
check_config
check_apex

OUT="${DATA_DIR}/figure3/sweep"
NUM_SAMPLES=10      # samples per prototype (50 protos × 10 = 500 sequences per run)
BATCH_SIZE=256

# Full τ and σ ranges used in the paper
TAU_VALUES=(0 50 100 150 200 300 500 700)
SIGMA_VALUES=(0.0 0.25 0.5 0.75 1.0)

echo "================================================================"
echo "Figure 3 — analog parameter sweep"
echo "Model:     ${CHECKPOINT}"
echo "Positives: ${POSITIVE_PROTOTYPES}"
echo "Negatives: ${NEGATIVE_PROTOTYPES}"
echo "Output:    ${OUT}"
printf "τ values:  %s\n" "${TAU_VALUES[*]}"
printf "σ values:  %s\n" "${SIGMA_VALUES[*]}"
echo "Combinations: $(( ${#TAU_VALUES[@]} * ${#SIGMA_VALUES[@]} )) per prototype set"
echo "================================================================"

# Run generation for one prototype set
run_sweep() {
    local proto_file="$1"
    local pt="$2"       # "positive" or "negative"

    local seq_dir="${OUT}/${pt}/sequences"
    local pred_dir="${OUT}/${pt}/predictions"
    mkdir -p "${seq_dir}" "${pred_dir}"

    local n_protos
    n_protos=$(grep -c "^>" "${proto_file}" || true)
    echo ""
    echo "━━━━ ${pt} prototypes (${n_protos} sequences) ━━━━"

    local total=$(( ${#TAU_VALUES[@]} * ${#SIGMA_VALUES[@]} ))
    local run=0

    for TAU in "${TAU_VALUES[@]}"; do
        for SIGMA in "${SIGMA_VALUES[@]}"; do
            (( run++ )) || true
            local BASE="${pt}_tau${TAU}_sigma${SIGMA}"
            local FASTA="${seq_dir}/${BASE}.fasta"
            local PRED="${pred_dir}/${BASE}-min.tsv"

            echo ""
            echo "  [${run}/${total}] τ=${TAU}, σ=${SIGMA}"

            # Generation
            if [ -f "${FASTA}" ] && [ -s "${FASTA}" ]; then
                echo "    sequences: skip (exists)"
            else
                cd "${OMEGAMP_DIR}"
                CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" AnalogConditional \
                    --checkpoint_path "${CHECKPOINT}" \
                    --analog_file "${proto_file}" \
                    --initial_timestep "${TAU}" \
                    --conditioning_closeness "${SIGMA}" \
                    --num_samples "${NUM_SAMPLES}" \
                    --batch_size "${BATCH_SIZE}" \
                    --output_fasta "${FASTA}" \
                    --conditioning_output_path "${seq_dir}/${BASE}_conditioning.pt"
                echo "    sequences: done"
            fi

            # APEX predictions
            if [ -f "${PRED}" ] && [ -s "${PRED}" ]; then
                echo "    predictions: skip (exists)"
            elif [ -f "${FASTA}" ] && [ -s "${FASTA}" ]; then
                run_apex "${FASTA}" "${PRED}"
                echo "    predictions: done"
            else
                echo "    predictions: skip (no FASTA)"
            fi
        done
    done
}

# Verify prototype files
for F in "${POSITIVE_PROTOTYPES}" "${NEGATIVE_PROTOTYPES}"; do
    if [ ! -f "${F}" ]; then
        echo "ERROR: prototype file not found: ${F}" >&2
        echo "  These files live in OMEGAMP_DIR/data/activity-data/." >&2
        echo "  Ensure OMEGAMP_DIR is set correctly: ${OMEGAMP_DIR}" >&2
        exit 1
    fi
done

run_sweep "${POSITIVE_PROTOTYPES}" "positive"
run_sweep "${NEGATIVE_PROTOTYPES}" "negative"

echo ""
echo "================================================================"
echo "Figure 3 sweep complete."
seq_count=$(find "${OUT}" -name "*.fasta" | wc -l)
pred_count=$(find "${OUT}" -name "*-min.tsv" | wc -l)
echo "  Sequence files:   ${seq_count}"
echo "  Prediction files: ${pred_count}"
echo "  Output: ${OUT}"
echo "================================================================"
