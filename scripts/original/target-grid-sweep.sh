#!/bin/bash

# ==============================================================================
# Generates samples across multiple 2D property grids:
#   1. Charge × Hydrophobicity (original)
#   2. Charge × Length
#   3. Hydrophobicity × Length
# ==============================================================================

set -e  # Exit on error

# Configuration
BASE_OUTPUT="experiments/results/files/phenotype-guided-pilot"
CHECKPOINT="models/generative_model.ckpt"
NUM_SAMPLES=250
BATCH_SIZE=250
CHARGE_GRID=(0 2 4 6 8 10 12)
HYDROPHOB_GRID=("-0.75" "-0.50" "-0.25" "0.0" "0.25" "0.5" "0.75")
LENGTH_GRID=(5 10 15 20 25 30 35)

# Create all necessary directories upfront
echo "Creating directory structure..."
mkdir -p "${BASE_OUTPUT}"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-charge-hydrophob"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-charge-length"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-hydrophob-length"

echo "=========================================="
echo "Checkpoint: ${CHECKPOINT}"
echo "Output: ${BASE_OUTPUT}"
echo "Samples per point: ${NUM_SAMPLES}"
echo "=========================================="

# ==============================================================================
# 1. CHARGE × HYDROPHOBICITY GRID (Original)
# ==============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Grid 1: Charge × Hydrophobicity"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"



for CHARGE in "${CHARGE_GRID[@]}"; do
    for HYDROPHOB in "${HYDROPHOB_GRID[@]}"; do
        HYDROPHOB_CLEAN=$(echo ${HYDROPHOB} | tr '.' 'p' | tr '-' 'm')
        echo "  Grid point: C=${CHARGE}, H=${HYDROPHOB}"

        OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-charge-hydrophob/c${CHARGE}_h${HYDROPHOB_CLEAN}"
        mkdir -p "${OUTPUT_DIR}"

        # Skip if already generated
        if [ -f "${OUTPUT_DIR}/samples.fasta" ] && [ -s "${OUTPUT_DIR}/samples.fasta" ]; then
            echo "    → Skipping (already exists)"
            continue
        fi

        CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
            "PartialConditional" \
            --checkpoint_path "${CHECKPOINT}" \
            --charge "${CHARGE}" \
            --length "-" \
            --hydrophobicity "${HYDROPHOB}" \
            --output_fasta "${OUTPUT_DIR}/samples.fasta" \
            --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE}
    done
done

# ==============================================================================
# 2. CHARGE × LENGTH GRID
# ==============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Grid 2: Charge × Length"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


for CHARGE in "${CHARGE_GRID[@]}"; do
    for LENGTH in "${LENGTH_GRID[@]}"; do
        echo "  Grid point: C=${CHARGE}, L=${LENGTH}"

        OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-charge-length/c${CHARGE}_l${LENGTH}"
        mkdir -p "${OUTPUT_DIR}"

        # Skip if already generated
        if [ -f "${OUTPUT_DIR}/samples.fasta" ] && [ -s "${OUTPUT_DIR}/samples.fasta" ]; then
            echo "    → Skipping (already exists)"
            continue
        fi

        CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
            "PartialConditional" \
            --checkpoint_path "${CHECKPOINT}" \
            --charge "${CHARGE}" \
            --length "${LENGTH}" \
            --hydrophobicity "-" \
            --output_fasta "${OUTPUT_DIR}/samples.fasta" \
            --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE}
    done
done

# ==============================================================================
# 3. HYDROPHOBICITY × LENGTH GRID
# ==============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Grid 3: Hydrophobicity × Length"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


for HYDROPHOB in "${HYDROPHOB_GRID[@]}"; do
    for LENGTH in "${LENGTH_GRID[@]}"; do
        HYDROPHOB_CLEAN=$(echo ${HYDROPHOB} | tr '.' 'p' | tr '-' 'm')
        echo "  Grid point: H=${HYDROPHOB}, L=${LENGTH}"

        OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-hydrophob-length/h${HYDROPHOB_CLEAN}_l${LENGTH}"
        mkdir -p "${OUTPUT_DIR}"

        # Skip if already generated
        if [ -f "${OUTPUT_DIR}/samples.fasta" ] && [ -s "${OUTPUT_DIR}/samples.fasta" ]; then
            echo "    → Skipping (already exists)"
            continue
        fi

        CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
            "PartialConditional" \
            --checkpoint_path "${CHECKPOINT}" \
            --charge "-" \
            --length "${LENGTH}" \
            --hydrophobicity "${HYDROPHOB}" \
            --output_fasta "${OUTPUT_DIR}/samples.fasta" \
            --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
            --num_samples ${NUM_SAMPLES} \
            --batch_size ${BATCH_SIZE}
    done
done


# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "=============================================="
echo "Grid Sweep Complete!"
echo "=============================================="

# Count generated files
count_grid1=$(find "${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-charge-hydrophob" -name "samples.fasta" 2>/dev/null | wc -l)
count_grid2=$(find "${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-charge-length" -name "samples.fasta" 2>/dev/null | wc -l)
count_grid3=$(find "${BASE_OUTPUT}/targeted-physicochemical/grid-sweep-hydrophob-length" -name "samples.fasta" 2>/dev/null | wc -l)

echo "Generated samples:"
echo "  • Charge × Hydrophobicity: ${count_grid1} grid points"
echo "  • Charge × Length:         ${count_grid2} grid points"
echo "  • Hydrophobicity × Length: ${count_grid3} grid points"
echo ""
echo "Total: $((count_grid1 + count_grid2 + count_grid3)) configurations"
echo "=============================================="