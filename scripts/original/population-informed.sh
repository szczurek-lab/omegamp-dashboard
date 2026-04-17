#!/bin/bash

# ==============================================================================
## PHENOTYPE-GUIDED DESIGN
# ==============================================================================


set -e  # Exit on error

BASE_OUTPUT="experiments/results/files/phenotype-guided-pilot"
CHECKPOINT="models/generative_model.ckpt"
GPU=2

# Create all necessary directories upfront
echo "Creating directory structure..."
mkdir -p "${BASE_OUTPUT}"
mkdir -p "${BASE_OUTPUT}/population-informed"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/single-property"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/multi-property"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/grid-sweep"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/impossible"
mkdir -p "${BASE_OUTPUT}/targeted-physicochemical/density-targets"

echo "=========================================="
echo "PILOT RUN: PHENOTYPE-GUIDED DESIGN"
echo "=========================================="

# ==============================================================================
# PART 1: POPULATION-INFORMED GENERATION
# ==============================================================================

echo ""
echo "=========================================="
echo "PART 1: POPULATION-INFORMED GENERATION"
echo "=========================================="

# --------------------------------------------------
# 1A: Generate large sample for density plots
# --------------------------------------------------
echo ""
echo "Generating unconditional samples..."
CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
    "Unconditional" \
    --checkpoint_path "${CHECKPOINT}" \
    --output_fasta "${BASE_OUTPUT}/population-informed/unconditional-samples.fasta" \
    --conditioning_output_path "${BASE_OUTPUT}/population-informed/unconditional-conditioning.pt" \
    --num_samples 5000 \
    --batch_size 1000

echo "Generating samples conditioned on HQ AMPs..."
CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
    "SubsetConditional" \
    --checkpoint_path "${CHECKPOINT}" \
    --subset_sequences "data/activity-data/curated-AMPs.fasta" \
    --output_fasta "${BASE_OUTPUT}/population-informed/subset-HQ-AMPs-samples.fasta" \
    --conditioning_output_path "${BASE_OUTPUT}/population-informed/subset-HQ-AMPs-conditioning.pt" \
    --num_samples 5000 \
    --batch_size 1000

# --------------------------------------------------
# 1C: Strain-specific conditioning
# --------------------------------------------------
echo ""
echo "Strain-specific conditioning..."

# Define strain-specific subsets (adjust paths to your actual data)
declare -A STRAIN_SUBSETS=(
    ["E_coli"]="data/activity-data/strain-species-data/species/escherichiacoli_positive.fasta"
    ["S_aureus"]="data/activity-data/strain-species-data/species/staphylococcusaureus_positive.fasta"
    ["P_aeruginosa"]="data/activity-data/strain-species-data/species/pseudomonasaeruginosa_positive.fasta"
    ["K_pneumoniae"]="data/activity-data/strain-species-data/species/klebsiellapneumoniae_positive.fasta"
    ["A_baummanii"]="data/activity-data/strain-species-data/species/acinetobacterbaumannii_positive.fasta"
)

for STRAIN in "${!STRAIN_SUBSETS[@]}"; do
    SUBSET_PATH="${STRAIN_SUBSETS[$STRAIN]}"

    if [ ! -f "${SUBSET_PATH}" ]; then
        echo "Warning: Subset file not found: ${SUBSET_PATH}. Skipping ${STRAIN}..."
        continue
    fi

    echo "Processing strain: ${STRAIN}"

    CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
        "SubsetConditional" \
        --checkpoint_path "${CHECKPOINT}" \
        --subset_sequences "${SUBSET_PATH}" \
        --output_fasta "${BASE_OUTPUT}/population-informed/strain-${STRAIN}-samples.fasta" \
        --conditioning_output_path "${BASE_OUTPUT}/population-informed/strain-${STRAIN}-conditioning.pt" \
        --num_samples 1000 \
        --batch_size 500
done

# ==============================================================================
# TARGETED PHYSICOCHEMICAL DESIGN
# ==============================================================================

echo ""
echo "=========================================="
echo "PART 2: TARGETED PHYSICOCHEMICAL DESIGN"
echo "=========================================="

# --------------------------------------------------
# 2A: Single Property Control
# --------------------------------------------------
echo ""
echo "Single property control experiments..."

# Charge sweep
CHARGE_VALUES=(0 2 4 6 8 10)
for CHARGE in "${CHARGE_VALUES[@]}"; do
    echo "  Testing charge = ${CHARGE}"
    OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/single-property/charge-${CHARGE}"
    mkdir -p "${OUTPUT_DIR}"

    CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
        "PartialConditional" \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "${CHARGE}" \
        --length "-" \
        --hydrophobicity "-" \
        --output_fasta "${OUTPUT_DIR}/samples.fasta" \
        --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
        --num_samples 500 \
        --batch_size 500
done

# Length sweep
LENGTH_VALUES=(10 15 20 25 30)
for LENGTH in "${LENGTH_VALUES[@]}"; do
    echo "  Testing length = ${LENGTH}"
    OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/single-property/length-${LENGTH}"
    mkdir -p "${OUTPUT_DIR}"

    CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
        "PartialConditional" \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "-" \
        --length "${LENGTH}" \
        --hydrophobicity "-" \
        --output_fasta "${OUTPUT_DIR}/samples.fasta" \
        --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
        --num_samples 500 \
        --batch_size 500
done

# Hydrophobicity sweep
HYDROPHOB_VALUES=("-0.5" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8")
for HYDROPHOB in "${HYDROPHOB_VALUES[@]}"; do
    echo "  Testing hydrophobicity = ${HYDROPHOB}"
    HYDROPHOB_CLEAN=$(echo ${HYDROPHOB} | tr '.' 'p' | tr '-' 'm')  # Clean for filename
    OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/single-property/hydrophob-${HYDROPHOB_CLEAN}"
    mkdir -p "${OUTPUT_DIR}"

    CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
        "PartialConditional" \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "-" \
        --length "-" \
        --hydrophobicity "${HYDROPHOB}" \
        --output_fasta "${OUTPUT_DIR}/samples.fasta" \
        --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
        --num_samples 500 \
        --batch_size 500
done

# --------------------------------------------------
# 2B: Multi-Property Control (for MAE calculation)
# --------------------------------------------------
echo ""
echo "2B: Multi-property control experiments..."

declare -a MULTI_TARGETS=(
    "2:10:-0.2:low_charge_short"
    "4:15:0.0:moderate_low"
    "6:20:0.3:moderate"
    "8:25:0.4:high_moderate"
    "10:30:0.6:high_long"
)

for TARGET in "${MULTI_TARGETS[@]}"; do
    IFS=':' read -r CHARGE LENGTH HYDROPHOB NAME <<< "$TARGET"
    echo "  Testing ${NAME}: C=${CHARGE}, L=${LENGTH}, H=${HYDROPHOB}"

    OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/multi-property/${NAME}"
    mkdir -p "${OUTPUT_DIR}"

    CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
        "PartialConditional" \
        --checkpoint_path "${CHECKPOINT}" \
        --charge "${CHARGE}" \
        --length "${LENGTH}" \
        --hydrophobicity "${HYDROPHOB}" \
        --output_fasta "${OUTPUT_DIR}/samples.fasta" \
        --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
        --num_samples 500 \
        --batch_size 500
done

# --------------------------------------------------
# Grid Sweep (for 2D heatmap)
# --------------------------------------------------
echo ""
echo "Grid sweep for heatmap..."

CHARGE_GRID=(0 2 4 6 8 10)
HYDROPHOB_GRID=("-0.5" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8")

for CHARGE in "${CHARGE_GRID[@]}"; do
    for HYDROPHOB in "${HYDROPHOB_GRID[@]}"; do
        HYDROPHOB_CLEAN=$(echo ${HYDROPHOB} | tr '.' 'p' | tr '-' 'm')
        echo "  Grid point: C=${CHARGE}, H=${HYDROPHOB}"

        OUTPUT_DIR="${BASE_OUTPUT}/targeted-physicochemical/grid-sweep/c${CHARGE}_h${HYDROPHOB_CLEAN}"
        mkdir -p "${OUTPUT_DIR}"

        CUDA_VISIBLE_DEVICES=${GPU} python project/scripts/inference/generate_samples.py \
            "PartialConditional" \
            --checkpoint_path "${CHECKPOINT}" \
            --charge "${CHARGE}" \
            --length "-" \
            --hydrophobicity "${HYDROPHOB}" \
            --output_fasta "${OUTPUT_DIR}/samples.fasta" \
            --conditioning_output_path "${OUTPUT_DIR}/conditioning.pt" \
            --num_samples 250 \
            --batch_size 250
    done
done
