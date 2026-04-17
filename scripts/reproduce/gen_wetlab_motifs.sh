#!/usr/bin/env bash
# =============================================================================
# Generate motif/template-conditioned peptides for wet-lab validation.
#
# Two passes:
#
#   Pass 1 — Template-only, LPS-binding motifs
#     Generates sequences with the GWKRKRFG octapeptide fixed at the N- or
#     C-terminus.  No analog prototype required; uses TemplateConditional with
#     AMP physicochemical constraints.
#
#     Templates:
#       LPS_Nterminal  GWKRKRFG_______   (motif at N-terminus, 15 aa)
#       LPS_Cterminal  _______GWKRKRFG   (motif at C-terminus, 15 aa)
#
#   Pass 2 — Template + analog conditioning, σ sweep
#     For five named prototype peptides each paired with a residue-level
#     template, runs AdvancedConditional(TemplateConditional,AnalogConditional,
#     PartialConditional) across σ ∈ {0, 0.25, 0.5, 0.75, 1.0} at τ=250.
#
#     Prototypes and matched templates:
#       cecropin    SWLSKTAKKLENSAKKRISEGIAIAIQGGPR
#                   ______________KK_____I___I_____
#       sarcotoxin  GWLKKIGKKIERVGQHTRDATIQGLGIAQQAANVAATAR
#                   _W_KK__________________________________
#       pa4         GFFALIPKIISSPLFKTLLSAVGSALSSSGGQE
#                   _______KII__P__K_LL_A____________
#       bZIP        MKQLEDKVEELLSKNYHLENEVARLKKLVG
#                   ___L___V__L___N__L___V__L___V_
#       LG21        LLPIVGNLLKSLLGWKRKRFG
#                   _____________GWKRKRFG
#
# Post-processing pipeline (each step applied when the script is present):
#
#   1. check_template.py
#      Filters the FASTA to sequences where all fixed template positions match.
#      CLI: python check_template.py <input.fasta> <TEMPLATE> <filtered.fasta>
#
#   2. synthesis_filtering.py
#      Evaluates synthesis feasibility and writes a CSV annotated with
#      is_peptide_probable (True/False).  Criteria:
#        • Hydrophilic abundance (KRHDESTNQY): 30 – 70 %
#        • Max consecutive Gly: ≤ 2
#        • Max run of any single bulky AA (FYWILVMRH): ≤ 1
#        • No known aggregation-prone motifs
#          (VLVL, KLLL, KLVFFA, NFGAIL, VVV, …  28 motifs total)
#        • Proline content: ≤ 20 %
#        • Cysteines: ≤ 1 (even-paired disulfides)
#        • Net charge (pH 7): 1.9 – 10.0
#      CLI: python synthesis_filtering.py --input_path <fasta> --output_path <csv>
#      Note: the CSV includes ALL sequences with the flag; filter downstream on
#      is_peptide_probable == True before ordering synthesis.
#
#   3. predict_sequences.py
#      Ranks the template-filtered FASTA with the full classifier ensemble
#      (broad AMP + species + strains + hemolytic), outputs top-500 by
#      ensemble vote.
#      CLI: python predict_sequences.py <fasta> --classifier all
#               --output_csv <csv> --get_top_k 500
#
#   4. calculate_similarity.py
#      Appends a "similarity" column (pairwise global alignment, %) to the
#      predictions CSV.  Modifies the CSV in place.
#      CLI: python calculate_similarity.py <ref_seq> <csv> [threshold]
#
# Writes:
#   data/wetlab/motif-lps/{LPS_Nterminal,LPS_Cterminal}-samples.fasta
#   data/wetlab/motif-lps/{LPS_Nterminal,LPS_Cterminal}-filtered.fasta
#   data/wetlab/motif-lps/{LPS_Nterminal,LPS_Cterminal}-synthesis.csv
#   data/wetlab/motif-lps/{LPS_Nterminal,LPS_Cterminal}-predictions.csv
#   data/wetlab/template-analog/{NAME}/sigma-{S}-samples.fasta
#   data/wetlab/template-analog/{NAME}/sigma-{S}-filtered.fasta
#   data/wetlab/template-analog/{NAME}/sigma-{S}-synthesis.csv
#   data/wetlab/template-analog/{NAME}/sigma-{S}-predictions.csv  (with similarity col)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/config.sh"
check_config

OUT="${DATA_DIR}/wetlab"
NUM_SAMPLES=50000
BATCH_SIZE=1000

# τ fixed for template+analog pass (high divergence, template enforces key residues)
TAU=250

# σ sweep for template+analog conditioning
SIGMA_VALUES=(0 0.25 0.5 0.75 1.0)

# AMP physicochemical constraints
LENGTH_RANGE="5:30"
CHARGE_RANGE="2:10"
HYDROPHOBICITY_RANGE="-0.5:0.8"

# Post-processing scripts
POSTPROC="${OMEGAMP_DIR}/project/scripts/postprocessing"
INFER="${OMEGAMP_DIR}/project/scripts/inference"

CHECK_TEMPLATE_SCRIPT="${POSTPROC}/check_template.py"
SYNTHESIS_SCRIPT="${POSTPROC}/synthesis_filtering.py"
PREDICT_SCRIPT="${INFER}/predict_sequences.py"
SIMILARITY_SCRIPT="${POSTPROC}/calculate_similarity.py"

echo "================================================================"
echo "Wet-lab motif/template generation"
echo "Model:  ${CHECKPOINT}"
echo "Output: ${OUT}"
echo "================================================================"

cd "${OMEGAMP_DIR}"

# =============================================================================
# Post-processing helpers
# =============================================================================

# Filter FASTA to sequences matching all fixed template positions.
# Positional args: input_fasta  template_string  output_fasta
# Skipped silently when check_template.py is not present.
maybe_filter() {
    local fasta="$1" template="$2" filtered="$3"
    if [ ! -f "${CHECK_TEMPLATE_SCRIPT}" ]; then return; fi
    if [ ! -f "${fasta}" ] || [ ! -s "${fasta}" ]; then return; fi
    python "${CHECK_TEMPLATE_SCRIPT}" "${fasta}" "${template}" "${filtered}"
    local n; n=$(grep -c '^>' "${filtered}" 2>/dev/null || echo 0)
    echo "    template filter: ${n} sequences match"
}

# Run synthesis feasibility annotation.
# Output CSV has one row per input sequence + is_peptide_probable column.
# Synthesis criteria (all must pass):
#   hydrophilic abundance 30-70%, consecutive Gly ≤2, bulky-AA run ≤1,
#   no aggregation motifs, proline ≤20%, cysteines ≤1, charge 1.9-10.0
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

# Score with all classifiers and keep top-500.
# Runs on the template-filtered FASTA when available; otherwise the full FASTA.
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

# Append pairwise similarity to prototype (modifies the CSV in place).
# Positional args: ref_sequence  csv_file
maybe_similarity() {
    local seq="$1" csv="$2"
    if [ ! -f "${SIMILARITY_SCRIPT}" ]; then return; fi
    if [ ! -f "${csv}" ]; then return; fi
    python "${SIMILARITY_SCRIPT}" "${seq}" "${csv}"
    echo "    similarity: done"
}

# =============================================================================
# Pass 1: Template-only, LPS-binding motifs (TemplateConditional + PartialConditional)
# =============================================================================
echo ""
echo "━━━━ Pass 1: LPS-binding template-only generation ━━━━"
echo "  Motif: GWKRKRFG (LPS-binding octapeptide)"
echo "  N-terminal: GWKRKRFG_______ — fixes motif at positions 1-8"
echo "  C-terminal: _______GWKRKRFG — fixes motif at positions 8-15"

LPS_OUT="${OUT}/motif-lps"
mkdir -p "${LPS_OUT}"

declare -A LPS_TEMPLATES=(
    ["LPS_Nterminal"]="GWKRKRFG_______"
    ["LPS_Cterminal"]="_______GWKRKRFG"
)

for MOTIF_NAME in "${!LPS_TEMPLATES[@]}"; do
    TEMPLATE="${LPS_TEMPLATES[$MOTIF_NAME]}"
    FASTA="${LPS_OUT}/${MOTIF_NAME}-samples.fasta"
    FILTERED="${LPS_OUT}/${MOTIF_NAME}-filtered.fasta"
    SYNTHESIS_CSV="${LPS_OUT}/${MOTIF_NAME}-synthesis.csv"
    PRED_CSV="${LPS_OUT}/${MOTIF_NAME}-predictions.csv"

    echo ""
    echo "  ${MOTIF_NAME}  template=${TEMPLATE}"

    # Generation
    if [ -f "${FASTA}" ] && [ -s "${FASTA}" ]; then
        echo "    sequences: skip (exists)"
    else
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" AdvancedConditional \
            --checkpoint_path "${CHECKPOINT}" \
            --motif "${TEMPLATE}" \
            --length "${LENGTH_RANGE}" \
            --charge "${CHARGE_RANGE}" \
            --hydrophobicity "${HYDROPHOBICITY_RANGE}" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE} \
            --output_fasta "${FASTA}" \
            --conditioning_output_path "${LPS_OUT}/${MOTIF_NAME}-conditioning.pt" \
            --advanced_conditioning_modes "TemplateConditional,PartialConditional"
        echo "    sequences: done"
    fi

    # Post-processing
    maybe_filter "${FASTA}" "${TEMPLATE}" "${FILTERED}"
    # Synthesis annotation on the template-matched set; fallback to full FASTA
    if [ -f "${FILTERED}" ] && [ -s "${FILTERED}" ]; then
        maybe_synthesize "${FILTERED}" "${SYNTHESIS_CSV}"
        maybe_predict    "${FILTERED}" "${PRED_CSV}"
    else
        maybe_synthesize "${FASTA}" "${SYNTHESIS_CSV}"
        maybe_predict    "${FASTA}"   "${PRED_CSV}"
    fi
    # No similarity step here — template-only runs have no single prototype seq
done

# =============================================================================
# Pass 2: Template + analog conditioning, σ sweep (τ=250 fixed)
# =============================================================================
echo ""
echo "━━━━ Pass 2: Template+analog conditioning (τ=${TAU}, σ sweep) ━━━━"
echo "  Conditioning modes: TemplateConditional + AnalogConditional + PartialConditional"
echo ""
echo "  Prototype / template pairs:"
echo "    cecropin    SWLSKTAKKLENSAKKRISEGIAIAIQGGPR"
echo "                ______________KK_____I___I_____"
echo "    sarcotoxin  GWLKKIGKKIERVGQHTRDATIQGLGIAQQAANVAATAR"
echo "                _W_KK__________________________________"
echo "    pa4         GFFALIPKIISSPLFKTLLSAVGSALSSSGGQE"
echo "                _______KII__P__K_LL_A____________"
echo "    bZIP        MKQLEDKVEELLSKNYHLENEVARLKKLVG"
echo "                ___L___V__L___N__L___V__L___V_"
echo "    LG21        LLPIVGNLLKSLLGWKRKRFG"
echo "                _____________GWKRKRFG"

# Prototype sequences and their matched residue templates.
# '_' = free position; letter = fixed amino acid (must appear in generated sequence).
declare -A PROTO_SEQ=(
    ["cecropin"]="SWLSKTAKKLENSAKKRISEGIAIAIQGGPR"
    ["sarcotoxin"]="GWLKKIGKKIERVGQHTRDATIQGLGIAQQAANVAATAR"
    ["pa4"]="GFFALIPKIISSPLFKTLLSAVGSALSSSGGQE"
    ["bZIP"]="MKQLEDKVEELLSKNYHLENEVARLKKLVG"
    ["LG21"]="LLPIVGNLLKSLLGWKRKRFG"
)

declare -A PROTO_TEMPLATE=(
    ["cecropin"]="______________KK_____I___I_____"
    ["sarcotoxin"]="_W_KK__________________________________"
    ["pa4"]="_______KII__P__K_LL_A____________"
    ["bZIP"]="___L___V__L___N__L___V__L___V_"
    ["LG21"]="_____________GWKRKRFG"
)

TA_OUT="${OUT}/template-analog"

for NAME in "${!PROTO_SEQ[@]}"; do
    SEQ="${PROTO_SEQ[$NAME]}"
    TEMPLATE="${PROTO_TEMPLATE[$NAME]}"
    NAME_OUT="${TA_OUT}/${NAME}"
    mkdir -p "${NAME_OUT}"

    echo ""
    echo "  ${NAME}"
    echo "    seq:      ${SEQ}"
    echo "    template: ${TEMPLATE}"

    for SIGMA in "${SIGMA_VALUES[@]}"; do
        LABEL="sigma-${SIGMA}"
        FASTA="${NAME_OUT}/${LABEL}-samples.fasta"
        FILTERED="${NAME_OUT}/${LABEL}-filtered.fasta"
        SYNTHESIS_CSV="${NAME_OUT}/${LABEL}-synthesis.csv"
        PRED_CSV="${NAME_OUT}/${LABEL}-predictions.csv"

        echo ""
        echo "    σ=${SIGMA}"

        # Generation
        if [ -f "${FASTA}" ] && [ -s "${FASTA}" ]; then
            echo "      sequences: skip (exists)"
        else
            CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" AdvancedConditional \
                --checkpoint_path "${CHECKPOINT}" \
                --analog "${SEQ}" \
                --motif "${TEMPLATE}" \
                --initial_timestep "${TAU}" \
                --conditioning_closeness "${SIGMA}" \
                --length "${LENGTH_RANGE}" \
                --charge "${CHARGE_RANGE}" \
                --hydrophobicity "${HYDROPHOBICITY_RANGE}" \
                --num_samples ${NUM_SAMPLES} \
                --batch_size ${BATCH_SIZE} \
                --output_fasta "${FASTA}" \
                --conditioning_output_path "${NAME_OUT}/${LABEL}-conditioning.pt" \
                --advanced_conditioning_modes "TemplateConditional,AnalogConditional,PartialConditional"
            echo "      sequences: done"
        fi

        # Post-processing
        maybe_filter "${FASTA}" "${TEMPLATE}" "${FILTERED}"
        # Prefer template-filtered set for downstream steps
        local_fasta="${FASTA}"
        if [ -f "${FILTERED}" ] && [ -s "${FILTERED}" ]; then
            local_fasta="${FILTERED}"
        fi
        maybe_synthesize "${local_fasta}" "${SYNTHESIS_CSV}"
        maybe_predict    "${local_fasta}" "${PRED_CSV}"
        # Similarity to prototype added to predictions CSV in place
        maybe_similarity "${SEQ}" "${PRED_CSV}"
    done
done

echo ""
echo "================================================================"
echo "Wet-lab motif/template generation complete."
lps_count=$(find "${LPS_OUT}" -name "*-samples.fasta" 2>/dev/null | wc -l)
ta_count=$(find "${TA_OUT}" -name "*-samples.fasta" 2>/dev/null | wc -l)
echo "  LPS-template FASTAs:      ${lps_count}"
echo "  Template+analog FASTAs:   ${ta_count}"
echo "  Output: ${OUT}"
echo "================================================================"
