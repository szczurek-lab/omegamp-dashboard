#!/bin/bash

# ==============================================================================
# Generate APEX Predictions for Phenotype-Guided Generation Samples
# ==============================================================================
#
# This script generates APEX MIC predictions for all generated samples:
# - Population-informed (unconditional, strain-specific)
# - Targeted physicochemical samples
#
# APEX provides MIC predictions for multiple bacterial strains
# ==============================================================================

set -e  # Exit on error

# Get absolute paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_INPUT="/home/pszymczak/omegamp/experiments/results/files/phenotype-guided-pilot"

BASE_OUTPUT="/home/pszymczak/omegamp/experiments/results/files/phenotype-guided-pilot/apex-predictions"
APEX_DIR="/raid/battleamp_root/battleamp-benchmark/battleamp/models/apex"
APEX_SCRIPT="${APEX_DIR}/predict.py"

# Create output directory
mkdir -p "${BASE_OUTPUT}"

echo "=========================================="
echo "APEX PREDICTIONS GENERATION"
echo "=========================================="
echo ""
echo "APEX directory: ${APEX_DIR}"
echo "APEX script: ${APEX_SCRIPT}"
echo "Input directory: ${BASE_INPUT}"
echo "Output directory: ${BASE_OUTPUT}"
echo ""

# Check if APEX script exists
if [ ! -f "${APEX_SCRIPT}" ]; then
    echo "❌ ERROR: APEX script not found at ${APEX_SCRIPT}"
    echo "Please update the APEX_DIR path in this script"
    exit 1
fi

# Check if APEX directory has required files
if [ ! -f "${APEX_DIR}/aaindex1.csv" ]; then
    echo "❌ ERROR: aaindex1.csv not found in ${APEX_DIR}"
    echo "Please ensure APEX directory is complete"
    exit 1
fi

# ==============================================================================
# PART 1: Population-Informed Generation Samples
# ==============================================================================
#
#echo "=========================================="
#echo "PART 1: POPULATION-INFORMED SAMPLES"
#echo "=========================================="
#echo ""
#
## Define samples to predict
#declare -a POPULATION_SAMPLES=(
#    "unconditional"
#    "subset-HQ-AMPs"
#    "strain-A_baummanii"
#    "strain-E_coli"
#    "strain-S_aureus"
#    "strain-P_aeruginosa"
#    "strain-K_pneumoniae"
#)
#
#for SAMPLE in "${POPULATION_SAMPLES[@]}"; do
#    INPUT_FASTA="${BASE_INPUT}/population-informed/${SAMPLE}-samples.fasta"
#    OUTPUT_TSV="${BASE_OUTPUT}/population-informed/${SAMPLE}-apex-predictions.tsv"
#
#    if [ ! -f "${INPUT_FASTA}" ]; then
#        echo "⚠️  Skipping ${SAMPLE}: FASTA not found"
#        continue
#    fi
#
#    echo "Predicting: ${SAMPLE}"
#    echo "  Input: ${INPUT_FASTA}"
#    echo "  Output: ${OUTPUT_TSV}"
#
#    mkdir -p "$(dirname ${OUTPUT_TSV})"
#
#    # Change to APEX directory and run prediction
#    cd "${APEX_DIR}"
#    python predict.py \
#        "${INPUT_FASTA}" \
#        "${OUTPUT_TSV}" \
#        all
#    cd "${SCRIPT_DIR}"
#
#    echo "  ✓ Complete"
#    echo ""
#done

# ==============================================================================
# PART 2: Targeted Physicochemical - Single Property Samples (optional subset)
# ==============================================================================

echo "=========================================="
echo "PART 2: TARGETED PHYSICOCHEMICAL SAMPLES"
echo "=========================================="


echo "2A: Single property samples..."

# Charge samples
for CHARGE in 0 2 4 6 8 10; do
    INPUT_FASTA="${BASE_INPUT}/targeted-physicochemical/single-property/charge-${CHARGE}/samples.fasta"
    OUTPUT_TSV="${BASE_OUTPUT}/targeted-physicochemical/single-property/charge-${CHARGE}-apex-predictions.tsv"

    if [ -f "${INPUT_FASTA}" ]; then
        echo "Predicting: charge-${CHARGE}"
        mkdir -p "$(dirname ${OUTPUT_TSV})"

        # Change to APEX directory before running
        cd "${APEX_DIR}"
        python predict.py "${INPUT_FASTA}" "${OUTPUT_TSV}" all
        cd "${SCRIPT_DIR}"

        echo "  ✓ Complete"
    fi
done

# Length samples
for LENGTH in 10 15 20 25 30; do
    INPUT_FASTA="${BASE_INPUT}/targeted-physicochemical/single-property/length-${LENGTH}/samples.fasta"
    OUTPUT_TSV="${BASE_OUTPUT}/targeted-physicochemical/single-property/length-${LENGTH}-apex-predictions.tsv"

    if [ -f "${INPUT_FASTA}" ]; then
        echo "Predicting: length-${LENGTH}"
        mkdir -p "$(dirname ${OUTPUT_TSV})"

        # Change to APEX directory before running
        cd "${APEX_DIR}"
        python predict.py "${INPUT_FASTA}" "${OUTPUT_TSV}" all
        cd "${SCRIPT_DIR}"

        echo "  ✓ Complete"
    fi
done

# Hydrophobicity samples
HYDROPHOB_VALUES=("-0.5" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8")
for HYDROPHOB in "${HYDROPHOB_VALUES[@]}"; do
    HYDROPHOB_CLEAN=$(echo ${HYDROPHOB} | tr '.' 'p' | tr '-' 'm')
    INPUT_FASTA="${BASE_INPUT}/targeted-physicochemical/single-property/hydrophob-${HYDROPHOB_CLEAN}/samples.fasta"
    OUTPUT_TSV="${BASE_OUTPUT}/targeted-physicochemical/single-property/hydrophob-${HYDROPHOB_CLEAN}-apex-predictions.tsv"

    if [ -f "${INPUT_FASTA}" ]; then
        echo "Predicting: hydrophob-${HYDROPHOB}"
        mkdir -p "$(dirname ${OUTPUT_TSV})"

        # Change to APEX directory before running
        cd "${APEX_DIR}"
        python predict.py "${INPUT_FASTA}" "${OUTPUT_TSV}" all
        cd "${SCRIPT_DIR}"

        echo "  ✓ Complete"
    fi
done

## ==============================================================================
## SUMMARY
## ==============================================================================
#
#echo ""
#echo "=========================================="
#echo "APEX PREDICTIONS COMPLETE"
#echo "=========================================="
#echo ""
#
## Count generated prediction files
#PRED_COUNT=$(find "${BASE_OUTPUT}" -name "*-apex-predictions.tsv" | wc -l)
#
#echo "Generated ${PRED_COUNT} APEX prediction files"
#echo ""
#echo "Prediction files saved in: ${BASE_OUTPUT}"
#echo ""
#echo "APEX output columns include MIC predictions for:"
#echo "  - E. coli (multiple strains)"
#echo "  - S. aureus (including MRSA)"
#echo "  - P. aeruginosa (PAO1, PA14)"
#echo "  - K. pneumoniae"
#echo "  - A. baumannii"
#echo "  - And many other bacterial species"
#echo ""
#echo "Next steps:"
#echo "  1. Use these MIC predictions for correlation analysis"
#echo "  2. Compare generated samples' predicted MIC patterns"
#echo "  3. Analyze if strain-conditioned samples show expected specificity"
#echo ""
#echo "Example: Load predictions in notebook:"
#echo "  df = pd.read_csv('${BASE_OUTPUT}/population-informed/strain-E_coli-apex-predictions.tsv', sep='\\t')"


# ==============================================================================
# PART 3: Grid Sweep Samples (for heatmap/3D surface plot)
# ==============================================================================


echo "=========================================="
echo "PART 3: GRID SWEEP SAMPLES"
echo "Charge vs Hydrophobicity Grid"
echo "=========================================="
echo ""

# Define grid parameters
CHARGE_VALUES=(0 2 4 6 8 10)
HYDROPHOB_VALUES=("-0.5" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8")

GRID_COUNT=0
TOTAL_GRID_POINTS=$((${#CHARGE_VALUES[@]} * ${#HYDROPHOB_VALUES[@]}))
for CHARGE in "${CHARGE_VALUES[@]}"; do
    for HYDROPHOB in "${HYDROPHOB_VALUES[@]}"; do
        ((GRID_COUNT++)) || true

        # Format hydrophobicity for filename (replace . with p, - with m)
        HYDROPHOB_CLEAN=$(echo ${HYDROPHOB} | tr '.' 'p' | tr '-' 'm')

        # Use the actual directory naming format: c{charge}_h{hydrophob}
        DIR_NAME="c${CHARGE}_h${HYDROPHOB_CLEAN}"

        INPUT_FASTA="${BASE_INPUT}/targeted-physicochemical/grid-sweep/${DIR_NAME}/samples.fasta"
        OUTPUT_TSV="${BASE_OUTPUT}/targeted-physicochemical/grid-sweep/${DIR_NAME}-apex-predictions.tsv"

        if [ ! -f "${INPUT_FASTA}" ]; then
            echo "⚠️  [${GRID_COUNT}/${TOTAL_GRID_POINTS}] Skipping charge=${CHARGE}, hydrophob=${HYDROPHOB}: FASTA not found"
            continue
        fi

        echo "[${GRID_COUNT}/${TOTAL_GRID_POINTS}] Predicting: charge=${CHARGE}, hydrophob=${HYDROPHOB}"
        echo "  Input: ${INPUT_FASTA}"
        echo "  Output: ${OUTPUT_TSV}"

        mkdir -p "$(dirname ${OUTPUT_TSV})"

        # Run APEX with full path and better error handling
        echo "  Running APEX..."
        cd "${APEX_DIR}"

        # Capture output and errors
        if python predict.py "${INPUT_FASTA}" "${OUTPUT_TSV}" all 2>&1; then
            cd "${SCRIPT_DIR}"
            echo "  ✓ Complete"
        else
            EXIT_CODE=$?
            cd "${SCRIPT_DIR}"
            echo "  ✗ APEX failed with exit code ${EXIT_CODE}"
            echo "  Check if APEX is working correctly"
            # Don't exit, continue with next sample
        fi
        echo ""
    done
done