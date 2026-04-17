#!/usr/bin/env bash
# =============================================================================
# Generate analog candidates for wet-lab validation.
#
# Two passes per target peptide:
#   1. Analog-only (τ=100, σ varied) — broad exploration across closeness values
#   2. Analog + property (τ=250, σ=0.25) — constrained to AMP physicochemical space
#
# Target peptides are the experimentally selected negative-control / scaffold
# sequences used in the wet-lab screen.
#
# Post-processing pipeline applied after each generation run:
#
#   1. synthesis_filtering.py
#      Annotates sequences with is_peptide_probable (True/False).  Criteria:
#        • Hydrophilic abundance (KRHDESTNQY): 30 – 70 %
#        • Max consecutive Gly: ≤ 2
#        • Max run of any single bulky AA (FYWILVMRH): ≤ 1
#        • No known aggregation-prone motifs (VLVL, KLLL, KLVFFA, …  28 total)
#        • Proline content: ≤ 20 %
#        • Cysteines: ≤ 1 (even-paired disulfides only)
#        • Net charge (pH 7): 1.9 – 10.0
#      Output CSV includes all sequences; filter on is_peptide_probable == True
#      before ordering synthesis.
#
#   2. predict_sequences.py
#      Runs the full classifier ensemble (broad AMP + species + strains +
#      hemolytic) and outputs the top-500 ranked sequences.
#
#   3. calculate_similarity.py  (analog+property pass only)
#      Appends pairwise similarity (%) to the prototype sequence.
#      Modifies the predictions CSV in place.
#
# Writes:
#   data/wetlab/analog-only/closeness-{S}/{NAME}-samples.fasta
#   data/wetlab/analog-only/closeness-{S}/{NAME}-synthesis.csv
#   data/wetlab/analog-only/closeness-{S}/{NAME}-predictions.csv
#   data/wetlab/analog-property/{NAME}-samples.fasta
#   data/wetlab/analog-property/{NAME}-synthesis.csv
#   data/wetlab/analog-property/{NAME}-predictions.csv  (with similarity col)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/config.sh"
check_config

OUT="${DATA_DIR}/wetlab"
NUM_SAMPLES=50000
BATCH_SIZE=1000

# Target peptides (sequences selected for wet-lab validation)
declare -A TARGETS=(
    ["GQ20"]="GQLNKFIKKAQRKFHEKFAK"
    ["Bacteroidin-2"]="MGNLVAIVGRPNVGKSTLFNRFH"
    ["Mammutin-1"]="WMTIHALKLSLSFKL"
    ["DBAASPS_20955"]="IGKLFKRIVERIKRFLRVLLRILR"
    ["BoCo1"]="NKIKFINKYVKKVQLKKILVKS"
    ["DBAASPS_20015"]="LGLFKLLLRLILKGFKL"
    ["DeNo1047"]="ALPSIIKGLLKKL"
)

# σ values for the analog-only closeness sweep
CLOSENESS_VALUES=(0 0.25 0.5 0.75 1)

# AMP physicochemical constraints for the analog+property pass
LENGTH_RANGE="5:30"
CHARGE_RANGE="2:10"
HYDROPHOBICITY_RANGE="-0.5:0.8"

# Post-processing scripts (all optional — steps are skipped when not present)
POSTPROC="${OMEGAMP_DIR}/project/scripts/postprocessing"
INFER="${OMEGAMP_DIR}/project/scripts/inference"

SYNTHESIS_SCRIPT="${POSTPROC}/synthesis_filtering.py"
PREDICT_SCRIPT="${INFER}/predict_sequences.py"
SIMILARITY_SCRIPT="${POSTPROC}/calculate_similarity.py"

echo "================================================================"
echo "Wet-lab analog generation"
echo "Model:   ${CHECKPOINT}"
echo "Targets: ${!TARGETS[*]}"
echo "Output:  ${OUT}"
echo "================================================================"

cd "${OMEGAMP_DIR}"

# =============================================================================
# Post-processing helpers
# =============================================================================

# Annotate FASTA with synthesis feasibility criteria → CSV with is_peptide_probable.
maybe_synthesize() {
    local fasta="$1" csv="$2"
    if [ ! -f "${SYNTHESIS_SCRIPT}" ]; then return; fi
    if [ ! -f "${fasta}" ] || [ ! -s "${fasta}" ]; then return; fi
    python "${SYNTHESIS_SCRIPT}" --input_path "${fasta}" --output_path "${csv}"
    local n_prob; n_prob=$(python -c "
import pandas as pd
df = pd.read_csv('${csv}')
print(df['is_peptide_probable'].sum())
" 2>/dev/null || echo "?")
    echo "    synthesis filter: ${n_prob} / $(grep -c '^>' "${fasta}") synthesizable"
}

# Score with all classifiers, keep top-500.
maybe_predict() {
    local fasta="$1" csv="$2"
    if [ ! -f "${PREDICT_SCRIPT}" ]; then return; fi
    if [ ! -f "${fasta}" ] || [ ! -s "${fasta}" ]; then return; fi
    CUDA_VISIBLE_DEVICES=${GPU} python "${PREDICT_SCRIPT}" \
        "${fasta}" \
        --classifier "all" \
        --output_csv "${csv}" \
        --get_top_k 500
    echo "    predictions: done (top-500)"
}

# Append pairwise similarity (%) to the prototype; modifies CSV in place.
maybe_similarity() {
    local seq="$1" csv="$2"
    if [ ! -f "${SIMILARITY_SCRIPT}" ]; then return; fi
    if [ ! -f "${csv}" ]; then return; fi
    python "${SIMILARITY_SCRIPT}" "${seq}" "${csv}"
    echo "    similarity: done"
}

# =============================================================================
# Pass 1: Analog-only, varying σ (τ=100 fixed)
# =============================================================================
echo ""
echo "━━━━ Pass 1: Analog-only sweep (τ=100, σ ∈ {${CLOSENESS_VALUES[*]}}) ━━━━"

for SIGMA in "${CLOSENESS_VALUES[@]}"; do
    OUT_DIR="${OUT}/analog-only/closeness-${SIGMA}"
    mkdir -p "${OUT_DIR}"

    for NAME in "${!TARGETS[@]}"; do
        SEQ="${TARGETS[$NAME]}"
        FASTA="${OUT_DIR}/${NAME}-samples.fasta"

        echo ""
        echo "  ${NAME}  σ=${SIGMA}"

        if [ -f "${FASTA}" ] && [ -s "${FASTA}" ]; then
            echo "    skip (exists)"
            continue
        fi

        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" AnalogConditional \
            --checkpoint_path "${CHECKPOINT}" \
            --analog "${SEQ}" \
            --initial_timestep 100 \
            --conditioning_closeness "${SIGMA}" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE} \
            --output_fasta "${FASTA}" \
            --conditioning_output_path "${OUT_DIR}/${NAME}-conditioning.pt" \
            --advanced_conditioning_modes "AnalogConditional"

        maybe_synthesize "${FASTA}" "${OUT_DIR}/${NAME}-synthesis.csv"
        maybe_predict    "${FASTA}" "${OUT_DIR}/${NAME}-predictions.csv"
    done
done

# =============================================================================
# Pass 2: Analog + physicochemical property constraints (τ=250, σ=0.25)
# =============================================================================
echo ""
echo "━━━━ Pass 2: Analog + property constraints (τ=250, σ=0.25) ━━━━"

OUT_DIR="${OUT}/analog-property"
mkdir -p "${OUT_DIR}"

for NAME in "${!TARGETS[@]}"; do
    SEQ="${TARGETS[$NAME]}"
    FASTA="${OUT_DIR}/${NAME}-samples.fasta"

    echo ""
    echo "  ${NAME}"

    if [ -f "${FASTA}" ] && [ -s "${FASTA}" ]; then
        echo "    skip (exists)"
        continue
    fi

    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" AdvancedConditional \
        --checkpoint_path "${CHECKPOINT}" \
        --analog "${SEQ}" \
        --initial_timestep 250 \
        --conditioning_closeness 0.25 \
        --length "${LENGTH_RANGE}" \
        --charge "${CHARGE_RANGE}" \
        --hydrophobicity "${HYDROPHOBICITY_RANGE}" \
        --num_samples ${NUM_SAMPLES} \
        --batch_size ${BATCH_SIZE} \
        --output_fasta "${FASTA}" \
        --conditioning_output_path "${OUT_DIR}/${NAME}-conditioning.pt" \
        --advanced_conditioning_modes "AnalogConditional,PartialConditional"

    maybe_synthesize "${FASTA}" "${OUT_DIR}/${NAME}-synthesis.csv"
    maybe_predict    "${FASTA}" "${OUT_DIR}/${NAME}-predictions.csv"
    # Similarity to prototype appended to predictions CSV in place
    maybe_similarity "${SEQ}" "${OUT_DIR}/${NAME}-predictions.csv"
done

echo ""
echo "================================================================"
echo "Wet-lab generation complete."
fasta_count=$(find "${OUT}" -name "*-samples.fasta" | wc -l)
echo "  FASTA files: ${fasta_count}"
echo "  Output: ${OUT}"
echo "================================================================"
