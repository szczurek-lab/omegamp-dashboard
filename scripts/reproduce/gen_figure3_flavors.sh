#!/usr/bin/env bash
# =============================================================================
# Generate the analog flavors comparison dataset for Figure 3.
#
# Six conditions per prototype set (positive / negative):
#   0. OmegAMP baselines — copied from data/figure2/benchmarks/
#   1. Prototypes         — the prototype FASTAs themselves
#   2. Analog only        — AnalogConditional
#   3. Analog + property  — AnalogConditional + PartialConditional
#   4. Analog + template  — AnalogConditional + TemplateConditional (motif)
#   5. Analog + template + property
#
# τ/σ parameters (chosen from sweep analysis):
#   Positive: τ=150, σ=1.0
#   Negative: τ=100, σ=1.0
#
# Writes:
#   data/figure3/flavors/{0_omegamp_baselines … 5_analog_template_property}/
#     sequences/  — generated FASTA files
#     predictions/  — APEX MIC TSVs
#     similarity/   — sequence similarity CSVs (if calculate_similarity.py exists)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/config.sh"
check_config
check_apex

OUT="${DATA_DIR}/figure3/flavors"
NUM_SAMPLES=10
BATCH_SIZE=256

POSITIVE_TAU=150;  POSITIVE_SIGMA=1.0
NEGATIVE_TAU=100;  NEGATIVE_SIGMA=1.0

LENGTH_RANGE="5:30"
CHARGE_RANGE="2:10"
HYDROPHOBICITY_RANGE="-0.5:0.8"

MOTIF="G____G"
MOTIF_SEED=42

# OmegAMP-U / OmegAMP-T baselines copied from the benchmarks folder
OMEGAMP_U_SRC="${DATA_DIR}/figure2/benchmarks/unconditional.fasta"
OMEGAMP_T_SRC="${DATA_DIR}/figure2/benchmarks/target-physicochemical.fasta"

# Optional: similarity script (skipped if not present)
SIMILARITY_SCRIPT="${OMEGAMP_DIR}/project/scripts/postprocessing/calculate_similarity.py"

echo "================================================================"
echo "Figure 3 — analog flavors comparison"
echo "Model:     ${CHECKPOINT}"
echo "Positives: ${POSITIVE_PROTOTYPES}  (τ=${POSITIVE_TAU}, σ=${POSITIVE_SIGMA})"
echo "Negatives: ${NEGATIVE_PROTOTYPES}  (τ=${NEGATIVE_TAU}, σ=${NEGATIVE_SIGMA})"
echo "Output:    ${OUT}"
echo "================================================================"

# Create all subdirectory trees upfront
for FLAVOR in 0_omegamp_baselines 1_prototypes 2_analog_only 3_analog_property 4_analog_template 5_analog_template_property; do
    mkdir -p "${OUT}/${FLAVOR}/predictions"
    mkdir -p "${OUT}/${FLAVOR}/similarity"
done

# ---------------------------------------------------------------------------
# Helper: run APEX on a FASTA if the TSV doesn't exist yet
# Usage: maybe_apex INPUT_FASTA OUTPUT_TSV
# ---------------------------------------------------------------------------
maybe_apex() {
    local fasta="$1" tsv="$2"
    if [ -f "${tsv}" ] && [ -s "${tsv}" ]; then
        echo "      predictions: skip (exists)"
    elif [ -f "${fasta}" ] && [ -s "${fasta}" ]; then
        run_apex "${fasta}" "${tsv}"
        echo "      predictions: done"
    else
        echo "      predictions: skip (no FASTA)"
    fi
}

# ---------------------------------------------------------------------------
# Helper: run similarity calculation if the script is available
# Usage: maybe_similarity PROTOTYPE_FASTA PREDICTION_CSV OUTPUT_CSV
# ---------------------------------------------------------------------------
maybe_similarity() {
    local proto="$1" pred_csv="$2" out_csv="$3"
    if [ -f "${SIMILARITY_SCRIPT}" ] && [ -f "${pred_csv}" ]; then
        python "${SIMILARITY_SCRIPT}" "${proto}" "${pred_csv}" --output "${out_csv}" 2>/dev/null || true
    fi
}

# ---------------------------------------------------------------------------
# 0. Baselines — copy from figure2/benchmarks
# ---------------------------------------------------------------------------
echo ""
echo "━━━━ Condition 0: OmegAMP baselines ━━━━"

for SRC_FILE in "${OMEGAMP_U_SRC}" "${OMEGAMP_T_SRC}"; do
    DEST_NAME=$(basename "${SRC_FILE}" .fasta)
    DEST="${OUT}/0_omegamp_baselines/${DEST_NAME}.fasta"
    if [ -f "${SRC_FILE}" ]; then
        cp "${SRC_FILE}" "${DEST}"
        echo "  copied: $(basename ${SRC_FILE})"
        maybe_apex "${DEST}" "${OUT}/0_omegamp_baselines/predictions/${DEST_NAME}-min.tsv"
    else
        echo "  WARNING: ${SRC_FILE} not found — run gen_figure2_benchmarks first" >&2
    fi
done

# ---------------------------------------------------------------------------
# Run conditions 1–5 for one prototype set
# ---------------------------------------------------------------------------
run_conditions() {
    local proto_file="$1"
    local pt="$2"
    local tau="$3"
    local sigma="$4"

    echo ""
    echo "======== ${pt} prototypes (τ=${tau}, σ=${sigma}) ========"

    cd "${OMEGAMP_DIR}"

    # ---- 1. Prototypes baseline ----
    echo ""
    echo "  Condition 1: prototypes"
    local DEST="${OUT}/1_prototypes/${pt}_prototypes.fasta"
    cp "${proto_file}" "${DEST}"
    maybe_apex "${DEST}" "${OUT}/1_prototypes/predictions/${pt}_prototypes-min.tsv"

    # ---- 2. Analog only ----
    echo ""
    echo "  Condition 2: analog only"
    local FASTA_2="${OUT}/2_analog_only/${pt}_analog_only.fasta"
    if [ ! -f "${FASTA_2}" ] || [ ! -s "${FASTA_2}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_ANALOG_SCRIPT}" AnalogConditional \
            --analog_file "${proto_file}" \
            --checkpoint_path "${CHECKPOINT}" \
            --initial_timestep "${tau}" \
            --conditioning_closeness "${sigma}" \
            --num_samples "${NUM_SAMPLES}" \
            --batch_size "${BATCH_SIZE}" \
            --output_fasta "${FASTA_2}" \
            --conditioning_output_path "${OUT}/2_analog_only/${pt}_conditioning.pt"
    else
        echo "    sequences: skip (exists)"
    fi
    maybe_apex "${FASTA_2}" "${OUT}/2_analog_only/predictions/${pt}_analog_only-min.tsv"

    # ---- 3. Analog + property ----
    echo ""
    echo "  Condition 3: analog + property"
    local FASTA_3="${OUT}/3_analog_property/${pt}_analog_property.fasta"
    if [ ! -f "${FASTA_3}" ] || [ ! -s "${FASTA_3}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_ANALOG_SCRIPT}" AdvancedConditional \
            --analog_file "${proto_file}" \
            --checkpoint_path "${CHECKPOINT}" \
            --length "${LENGTH_RANGE}" \
            --charge "${CHARGE_RANGE}" \
            --hydrophobicity "${HYDROPHOBICITY_RANGE}" \
            --initial_timestep "${tau}" \
            --conditioning_closeness "${sigma}" \
            --num_samples "${NUM_SAMPLES}" \
            --batch_size "${BATCH_SIZE}" \
            --output_fasta "${FASTA_3}" \
            --conditioning_output_path "${OUT}/3_analog_property/${pt}_conditioning.pt" \
            --advanced_conditioning_modes "AnalogConditional,PartialConditional"
    else
        echo "    sequences: skip (exists)"
    fi
    maybe_apex "${FASTA_3}" "${OUT}/3_analog_property/predictions/${pt}_analog_property-min.tsv"

    # ---- 4. Analog + template (motif) ----
    echo ""
    echo "  Condition 4: analog + template"
    local FASTA_4="${OUT}/4_analog_template/${pt}_analog_template.fasta"
    if [ ! -f "${FASTA_4}" ] || [ ! -s "${FASTA_4}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_ANALOG_SCRIPT}" AdvancedConditional \
            --analog_file "${proto_file}" \
            --checkpoint_path "${CHECKPOINT}" \
            --motif "${MOTIF}" \
            --motif_seed "${MOTIF_SEED}" \
            --initial_timestep "${tau}" \
            --conditioning_closeness "${sigma}" \
            --num_samples "${NUM_SAMPLES}" \
            --batch_size "${BATCH_SIZE}" \
            --output_fasta "${FASTA_4}" \
            --conditioning_output_path "${OUT}/4_analog_template/${pt}_conditioning.pt" \
            --advanced_conditioning_modes "TemplateConditional,AnalogConditional"
    else
        echo "    sequences: skip (exists)"
    fi
    maybe_apex "${FASTA_4}" "${OUT}/4_analog_template/predictions/${pt}_analog_template-min.tsv"

    # ---- 5. Analog + template + property ----
    echo ""
    echo "  Condition 5: analog + template + property"
    local FASTA_5="${OUT}/5_analog_template_property/${pt}_analog_template_property.fasta"
    if [ ! -f "${FASTA_5}" ] || [ ! -s "${FASTA_5}" ]; then
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_ANALOG_SCRIPT}" AdvancedConditional \
            --analog_file "${proto_file}" \
            --checkpoint_path "${CHECKPOINT}" \
            --motif "${MOTIF}" \
            --motif_seed "${MOTIF_SEED}" \
            --length "${LENGTH_RANGE}" \
            --charge "${CHARGE_RANGE}" \
            --hydrophobicity "${HYDROPHOBICITY_RANGE}" \
            --initial_timestep "${tau}" \
            --conditioning_closeness "${sigma}" \
            --num_samples "${NUM_SAMPLES}" \
            --batch_size "${BATCH_SIZE}" \
            --output_fasta "${FASTA_5}" \
            --conditioning_output_path "${OUT}/5_analog_template_property/${pt}_conditioning.pt" \
            --advanced_conditioning_modes "TemplateConditional,AnalogConditional,PartialConditional"
    else
        echo "    sequences: skip (exists)"
    fi
    maybe_apex "${FASTA_5}" "${OUT}/5_analog_template_property/predictions/${pt}_analog_template_property-min.tsv"
}

# Verify prototype files
for F in "${POSITIVE_PROTOTYPES}" "${NEGATIVE_PROTOTYPES}"; do
    if [ ! -f "${F}" ]; then
        echo "ERROR: prototype file not found: ${F}" >&2
        exit 1
    fi
done

run_conditions "${POSITIVE_PROTOTYPES}" "positive" "${POSITIVE_TAU}" "${POSITIVE_SIGMA}"
run_conditions "${NEGATIVE_PROTOTYPES}" "negative" "${NEGATIVE_TAU}" "${NEGATIVE_SIGMA}"

echo ""
echo "================================================================"
echo "Figure 3 flavors comparison complete."
fasta_count=$(find "${OUT}" -name "*.fasta" | wc -l)
pred_count=$(find "${OUT}" -name "*-min.tsv" | wc -l)
echo "  FASTA files:      ${fasta_count}"
echo "  Prediction files: ${pred_count}"
echo "  Output: ${OUT}"
echo "================================================================"
