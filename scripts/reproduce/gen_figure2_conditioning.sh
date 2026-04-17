#!/usr/bin/env bash
# =============================================================================
# Generate all conditioning samples needed for Figure 2.
#
# Produces:
#   data/figure2/conditioning/single-property/{charge-N,length-N,hydrophob-X}/samples.fasta
#   data/figure2/conditioning/multi-property/{condition_name}/samples.fasta
#   data/figure2/conditioning/grid-sweep-charge-hydrophob/c{C}_h{H}/samples.fasta
#   data/figure2/conditioning/grid-sweep-charge-length/c{C}_l{L}/samples.fasta
#   data/figure2/conditioning/grid-sweep-hydrophob-length/h{H}_l{L}/samples.fasta
#
# Runtime: ~2–4 h on a single GPU depending on hardware.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/config.sh"
check_config

OUT="${DATA_DIR}/figure2/conditioning"
NUM_SAMPLES=250
BATCH_SIZE=250

# Grid axes (must match what the figure2 notebook expects)
CHARGE_GRID=(0 2 4 6 8 10 12)
HYDROPHOB_GRID=("-0.75" "-0.50" "-0.25" "0.0" "0.25" "0.5" "0.75")
LENGTH_GRID=(5 10 15 20 25 30 35)

# Single-property axes (panel B of figure 2)
CHARGE_SINGLE=(0 2 4 6 8 10)
HYDROPHOB_SINGLE=("-0.5" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8")
LENGTH_SINGLE=(10 15 20 25 30)

# Multi-property conditions (panel B)
declare -A MULTI_CHARGE=( [low_charge_short]=2  [moderate_low]=4  [moderate]=6  [high_moderate]=8  [high_long]=10 )
declare -A MULTI_LENGTH=( [low_charge_short]=10 [moderate_low]=15 [moderate]=20 [high_moderate]=25 [high_long]=30 )
declare -A MULTI_HYDRO=(  [low_charge_short]=-0.2 [moderate_low]=0.0 [moderate]=0.3 [high_moderate]=0.4 [high_long]=0.6 )

echo "================================================================"
echo "Figure 2 — conditioning data generation"
echo "Model:  ${CHECKPOINT}"
echo "Output: ${OUT}"
echo "Samples per point: ${NUM_SAMPLES}"
echo "================================================================"

# All generation runs from the OMEGAMP_DIR so Hydra can find its config
cd "${OMEGAMP_DIR}"

# Helper: skip grid point if output already exists
_skip() { [ -f "$1" ] && [ -s "$1" ] && echo "    → skip (exists)" && return 0 || return 1; }

# =============================================================================
# 1. Single-property conditioning
# =============================================================================
echo ""
echo "━━━━ 1. Single-property conditioning ━━━━"

echo "  Charge..."
for C in "${CHARGE_SINGLE[@]}"; do
    DIR="${OUT}/single-property/charge-${C}"
    mkdir -p "${DIR}"
    _skip "${DIR}/samples.fasta" && continue
    echo "    charge=${C}"
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" PartialConditional \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "${C}" \
        --length "-" \
        --hydrophobicity "-" \
        --output_fasta "${DIR}/samples.fasta" \
        --conditioning_output_path "${DIR}/conditioning.pt" \
        --num_samples ${NUM_SAMPLES} \
        --batch_size ${BATCH_SIZE}
done

echo "  Length..."
for L in "${LENGTH_SINGLE[@]}"; do
    DIR="${OUT}/single-property/length-${L}"
    mkdir -p "${DIR}"
    _skip "${DIR}/samples.fasta" && continue
    echo "    length=${L}"
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" PartialConditional \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "-" \
        --length "${L}" \
        --hydrophobicity "-" \
        --output_fasta "${DIR}/samples.fasta" \
        --conditioning_output_path "${DIR}/conditioning.pt" \
        --num_samples ${NUM_SAMPLES} \
        --batch_size ${BATCH_SIZE}
done

echo "  Hydrophobicity..."
for H in "${HYDROPHOB_SINGLE[@]}"; do
    HC=$(fmt_hydro "${H}")
    DIR="${OUT}/single-property/hydrophob-${HC}"
    mkdir -p "${DIR}"
    _skip "${DIR}/samples.fasta" && continue
    echo "    hydrophobicity=${H}"
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" PartialConditional \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "-" \
        --length "-" \
        --hydrophobicity "${H}" \
        --output_fasta "${DIR}/samples.fasta" \
        --conditioning_output_path "${DIR}/conditioning.pt" \
        --num_samples ${NUM_SAMPLES} \
        --batch_size ${BATCH_SIZE}
done

# =============================================================================
# 2. Multi-property conditioning
# =============================================================================
echo ""
echo "━━━━ 2. Multi-property conditioning ━━━━"

for COND_NAME in low_charge_short moderate_low moderate high_moderate high_long; do
    C="${MULTI_CHARGE[$COND_NAME]}"
    L="${MULTI_LENGTH[$COND_NAME]}"
    H="${MULTI_HYDRO[$COND_NAME]}"
    DIR="${OUT}/multi-property/${COND_NAME}"
    mkdir -p "${DIR}"
    _skip "${DIR}/samples.fasta" && { echo "    ${COND_NAME} → skip"; continue; }
    echo "  ${COND_NAME}: C=${C}, L=${L}, H=${H}"
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" PartialConditional \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "${C}" \
        --length "${L}" \
        --hydrophobicity "${H}" \
        --output_fasta "${DIR}/samples.fasta" \
        --conditioning_output_path "${DIR}/conditioning.pt" \
        --num_samples ${NUM_SAMPLES} \
        --batch_size ${BATCH_SIZE}
done

# =============================================================================
# 3. Grid sweep — Charge × Hydrophobicity
# =============================================================================
echo ""
echo "━━━━ 3. Grid sweep: Charge × Hydrophobicity ━━━━"

for C in "${CHARGE_GRID[@]}"; do
    for H in "${HYDROPHOB_GRID[@]}"; do
        HC=$(fmt_hydro "${H}")
        DIR="${OUT}/grid-sweep-charge-hydrophob/c${C}_h${HC}"
        mkdir -p "${DIR}"
        _skip "${DIR}/samples.fasta" && continue
        echo "  C=${C}, H=${H}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" PartialConditional \
            --checkpoint_path "${CHECKPOINT}" \
            --charge "${C}" \
            --length "-" \
            --hydrophobicity "${H}" \
            --output_fasta "${DIR}/samples.fasta" \
            --conditioning_output_path "${DIR}/conditioning.pt" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE}
    done
done

# =============================================================================
# 4. Grid sweep — Charge × Length
# =============================================================================
echo ""
echo "━━━━ 4. Grid sweep: Charge × Length ━━━━"

for C in "${CHARGE_GRID[@]}"; do
    for L in "${LENGTH_GRID[@]}"; do
        DIR="${OUT}/grid-sweep-charge-length/c${C}_l${L}"
        mkdir -p "${DIR}"
        _skip "${DIR}/samples.fasta" && continue
        echo "  C=${C}, L=${L}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" PartialConditional \
            --checkpoint_path "${CHECKPOINT}" \
            --charge "${C}" \
            --length "${L}" \
            --hydrophobicity "-" \
            --output_fasta "${DIR}/samples.fasta" \
            --conditioning_output_path "${DIR}/conditioning.pt" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE}
    done
done

# =============================================================================
# 5. Grid sweep — Hydrophobicity × Length
# =============================================================================
echo ""
echo "━━━━ 5. Grid sweep: Hydrophobicity × Length ━━━━"

for H in "${HYDROPHOB_GRID[@]}"; do
    for L in "${LENGTH_GRID[@]}"; do
        HC=$(fmt_hydro "${H}")
        DIR="${OUT}/grid-sweep-hydrophob-length/h${HC}_l${L}"
        mkdir -p "${DIR}"
        _skip "${DIR}/samples.fasta" && continue
        echo "  H=${H}, L=${L}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" PartialConditional \
            --checkpoint_path "${CHECKPOINT}" \
            --charge "-" \
            --length "${L}" \
            --hydrophobicity "${H}" \
            --output_fasta "${DIR}/samples.fasta" \
            --conditioning_output_path "${DIR}/conditioning.pt" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE}
    done
done

# =============================================================================
echo ""
echo "================================================================"
echo "Figure 2 conditioning generation complete."
fasta_count=$(find "${OUT}" -name "samples.fasta" | wc -l)
echo "Total FASTA files: ${fasta_count}"
echo "================================================================"
