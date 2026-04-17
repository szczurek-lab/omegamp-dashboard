#!/bin/bash

# Parameter Sweep for Analog Generation (τ and σ)
# This script runs analog generation across a grid of parameters for the Cell manuscript

set -e  # Exit on error

# ============================================================================
# CONFIGURATION
# ============================================================================

# Model and data paths
CHECKPOINT_PATH="models/generative_model.ckpt"

# Prototype files - START WITH SAMPLES for initial testing
POSITIVE_PROTOTYPES="data/activity-data/curated-AMPs_samples.fasta"
NEGATIVE_PROTOTYPES="data/activity-data/curated-Non-AMPs_samples.fasta"

# For full dataset, uncomment these:
# POSITIVE_PROTOTYPES="data/activity-data/curated-AMPs.fasta"
# NEGATIVE_PROTOTYPES="data/activity-data/curated-Non-AMPs.fasta"

# Output directory structure
OUTPUT_DIR="experiments/results/files/analog-parameter-sweep"
POSITIVE_DIR="${OUTPUT_DIR}/positive"
NEGATIVE_DIR="${OUTPUT_DIR}/negative"

# Create output directories
mkdir -p "${POSITIVE_DIR}/sequences"
mkdir -p "${POSITIVE_DIR}/conditioning"
mkdir -p "${POSITIVE_DIR}/metadata"
mkdir -p "${NEGATIVE_DIR}/sequences"
mkdir -p "${NEGATIVE_DIR}/conditioning"
mkdir -p "${NEGATIVE_DIR}/metadata"

# Sampling parameters
NUM_SAMPLES=10  # SAMPLES PER PROTOTYPE (not total)
                  # e.g., 50 prototypes × 5000 samples = 250,000 total sequences
BATCH_SIZE=256

# Parameter ranges for sweep
# τ (tau): exploration strength (initial_timestep) - higher means more exploration
#TAU_VALUES=(0 100 300 500 700)
TAU_VALUES=(50 150 200)

# σ (sigma): property relaxation (conditioning_closeness) - higher means more deviation
SIGMA_VALUES=(0.0 0.25 0.5 0.75 1.0)

# ============================================================================
# FUNCTIONS
# ============================================================================

run_sweep() {
    local prototype_file=$1
    local output_dir=$2
    local prototype_type=$3

    # Count prototypes in file
    local n_prototypes=$(grep -c "^>" "${prototype_file}" || echo "0")

    echo "========================================"
    echo "Running sweep for ${prototype_type} prototypes"
    echo "========================================"
    echo "Prototype file: ${prototype_file}"
    echo "Number of prototypes: ${n_prototypes}"
    echo "Samples per prototype: ${NUM_SAMPLES}"
    echo "Total samples per run: $((n_prototypes * NUM_SAMPLES))"
    echo ""

    local total_runs=$((${#TAU_VALUES[@]} * ${#SIGMA_VALUES[@]}))
    local current_run=0

    for tau in "${TAU_VALUES[@]}"; do
        for sigma in "${SIGMA_VALUES[@]}"; do
            current_run=$((current_run + 1))

            echo ""
            echo "----------------------------------------"
            echo "Run ${current_run}/${total_runs}: τ=${tau}, σ=${sigma}"
            echo "----------------------------------------"

            # Create filenames with parameters
            local base_name="${prototype_type}_tau${tau}_sigma${sigma}"
            local fasta_out="${output_dir}/sequences/${base_name}.fasta"
            local cond_out="${output_dir}/conditioning/${base_name}_conditioning.pt"
            local meta_out="${output_dir}/metadata/${base_name}_metadata.pt"

            # Run generation
            python project/scripts/inference/generate_samples.py AnalogConditional \
                --checkpoint_path "${CHECKPOINT_PATH}" \
                --analog_file "${prototype_file}" \
                --initial_timestep "${tau}" \
                --conditioning_closeness "${sigma}" \
                --num_samples "${NUM_SAMPLES}" \
                --batch_size "${BATCH_SIZE}" \
                --output_fasta "${fasta_out}" \
                --conditioning_output_path "${cond_out}" \
                --metadata_output_path "${meta_out}"

            echo "✓ Completed: ${base_name}"
            echo "  Generated: $((n_prototypes * NUM_SAMPLES)) sequences"
        done
    done

    echo ""
    echo "========================================"
    echo "Completed ${prototype_type} prototypes sweep"
    echo "========================================"
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

echo "================================================================"
echo "ANALOG GENERATION PARAMETER SWEEP"
echo "================================================================"
echo "Model: ${CHECKPOINT_PATH}"
echo "Positive prototypes: ${POSITIVE_PROTOTYPES}"
echo "Negative prototypes: ${NEGATIVE_PROTOTYPES}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Parameter ranges:"
echo "  τ (exploration): ${TAU_VALUES[@]}"
echo "  σ (relaxation):  ${SIGMA_VALUES[@]}"
echo ""
echo "Sampling configuration:"
echo "  Samples per prototype: ${NUM_SAMPLES}"
echo "  Batch size: ${BATCH_SIZE}"
echo ""
echo "Total parameter combinations: $((${#TAU_VALUES[@]} * ${#SIGMA_VALUES[@]}))"
echo "================================================================"
echo ""

# Check if prototype files exist
if [ ! -f "${POSITIVE_PROTOTYPES}" ]; then
    echo "ERROR: Positive prototype file not found: ${POSITIVE_PROTOTYPES}"
    echo ""
    echo "Please run: python scripts/analysis/create_prototype_samples.py"
    exit 1
fi

if [ ! -f "${NEGATIVE_PROTOTYPES}" ]; then
    echo "ERROR: Negative prototype file not found: ${NEGATIVE_PROTOTYPES}"
    echo ""
    echo "Please run: python scripts/analysis/create_prototype_samples.py"
    exit 1
fi

# Prompt for confirmation
read -p "Press Enter to start the sweep, or Ctrl+C to cancel..."

# Record start time
start_time=$(date +%s)

# Run sweeps
run_sweep "${POSITIVE_PROTOTYPES}" "${POSITIVE_DIR}" "positive"
run_sweep "${NEGATIVE_PROTOTYPES}" "${NEGATIVE_DIR}" "negative"

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo ""
echo "================================================================"
echo "SWEEP COMPLETED"
echo "================================================================"
echo "Total time: ${hours}h ${minutes}m ${seconds}s"
echo "Results saved to: ${OUTPUT_DIR}"
echo "================================================================"

# Save sweep configuration for reference
cat > "${OUTPUT_DIR}/sweep_config.txt" << EOF
Parameter Sweep Configuration
==============================
Date: $(date)
Model: ${CHECKPOINT_PATH}
Positive prototypes: ${POSITIVE_PROTOTYPES}
Negative prototypes: ${NEGATIVE_PROTOTYPES}

Parameters:
-----------
τ (tau) values: ${TAU_VALUES[@]}
σ (sigma) values: ${SIGMA_VALUES[@]}

Sampling:
---------
Samples per prototype: ${NUM_SAMPLES}
Total samples per run: varies by number of prototypes in file
Batch size: ${BATCH_SIZE}

Duration: ${hours}h ${minutes}m ${seconds}s
EOF

echo ""
echo "Configuration saved to: ${OUTPUT_DIR}/sweep_config.txt"