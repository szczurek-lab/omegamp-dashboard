#!/bin/bash

# ============================================================================
# ANALOG GENERATION FLAVORS COMPARISON
# ============================================================================
# Uses different τ/σ settings for positive vs negative prototypes:
#   - Positive: τ=50, σ=0.25  (stay close to known AMPs)
#   - Negative: τ=100, σ=1.0 (need flexibility to escape non-AMP space)
#
# Conditions evaluated:
#   0. OmegAMP baselines (copied from benchmark results)
#      - OmegAMP-U: Unconditional generation
#      - OmegAMP-T: Targeted physicochemical generation
#   1. Prototypes (baseline)
#   2. Analog only
#   3. Analog + Targeted physicochemical
#   4. Analog + Template (motif insertion)
#   5. Analog + Template + Targeted physicochemical
#
# Template conditioning uses dynamic motif insertion: a conserved motif
# pattern (e.g., "G____G") is inserted at random positions within each
# prototype, generating length-matched templates automatically.
# ============================================================================

set -e  # Exit on error

# ============================================================================
# CONFIGURATION
# ============================================================================

CHECKPOINT_PATH="models/generative_model.ckpt"

# Prototype files
POSITIVE_PROTOTYPES="data/activity-data/curated-AMPs_samples.fasta"
NEGATIVE_PROTOTYPES="data/activity-data/curated-Non-AMPs_samples.fasta"

# Baseline files (from benchmark results)
BENCHMARK_BASE="results/benchmark-generative-models"
OMEGAMP_U_SOURCE="${BENCHMARK_BASE}/unconditional.fasta"
OMEGAMP_T_SOURCE="${BENCHMARK_BASE}/target-physicochemical.fasta"

# Output directory
OUTPUT_DIR="experiments/results/files/analog-flavors-comparison"

# Sampling parameters
NUM_SAMPLES=10        # Samples per prototype
BATCH_SIZE=256

# Optimal parameters (from parameter sweep analysis)
# Positive: prioritize staying close to known AMPs
POSITIVE_TAU=150
POSITIVE_SIGMA=1.0

# Negative: need more flexibility to find functional sequences
NEGATIVE_TAU=100
NEGATIVE_SIGMA=1.0

# Property conditioning ranges (for targeted physicochemical)
LENGTH_RANGE="5:30"
CHARGE_RANGE="2:10"
HYDROPHOBICITY_RANGE="-0.5:0.8"

# Template motif for dynamic insertion
# Glycine spacer motif - will be inserted at random positions in each prototype
MOTIF="G____G"
MOTIF_SEED=42  # For reproducibility

# Path to generation scripts
GENERATE_SCRIPT="project/scripts/inference/generate_analog_samples.py"

# ============================================================================
# FUNCTIONS
# ============================================================================

create_directories() {
    local base_dir=$1
    mkdir -p "${base_dir}/0_omegamp_baselines"
    mkdir -p "${base_dir}/1_prototypes"
    mkdir -p "${base_dir}/2_analog_only"
    mkdir -p "${base_dir}/3_analog_property"
    mkdir -p "${base_dir}/4_analog_template"
    mkdir -p "${base_dir}/5_analog_template_property"
}



# Condition 0: Copy OmegAMP baselines
copy_omegamp_baselines() {
    local output_dir=$1

    echo "----------------------------------------"
    echo "Condition 0: Copying OmegAMP baselines"
    echo "----------------------------------------"

    # Copy OmegAMP-U (Unconditional)
    if [ -f "${OMEGAMP_U_SOURCE}" ]; then
        cp "${OMEGAMP_U_SOURCE}" "${output_dir}/0_omegamp_baselines/OmegAMP-U.fasta"
        echo "✓ Copied OmegAMP-U (unconditional) from ${OMEGAMP_U_SOURCE}"
    else
        echo "⚠ WARNING: OmegAMP-U source not found: ${OMEGAMP_U_SOURCE}"
    fi

    # Copy OmegAMP-T (Targeted physicochemical)
    if [ -f "${OMEGAMP_T_SOURCE}" ]; then
        cp "${OMEGAMP_T_SOURCE}" "${output_dir}/0_omegamp_baselines/OmegAMP-T.fasta"
        echo "✓ Copied OmegAMP-T (targeted) from ${OMEGAMP_T_SOURCE}"
    else
        echo "⚠ WARNING: OmegAMP-T source not found: ${OMEGAMP_T_SOURCE}"
    fi

    echo "✓ OmegAMP baselines copied to ${output_dir}/0_omegamp_baselines/"
}

# Copy prototypes as baseline (Condition 1: Prototypes)
copy_prototypes() {
    local prototype_file=$1
    local output_dir=$2
    local proto_type=$3

    echo "----------------------------------------"
    echo "Condition 1: Copying ${proto_type} prototypes as baseline"
    echo "----------------------------------------"

    cp "${prototype_file}" "${output_dir}/1_prototypes/${proto_type}_prototypes.fasta"
    echo "✓ Copied prototypes to ${output_dir}/1_prototypes/"
}

# Condition 2: Analog only
run_analog_only() {
    local prototype_file=$1
    local output_dir=$2
    local proto_type=$3
    local tau=$4
    local sigma=$5

    echo "----------------------------------------"
    echo "Condition 2: Analog only (${proto_type}, τ=${tau}, σ=${sigma})"
    echo "----------------------------------------"

    python "${GENERATE_SCRIPT}" AnalogConditional \
        --analog_file "${prototype_file}" \
        --checkpoint_path "${CHECKPOINT_PATH}" \
        --initial_timestep "${tau}" \
        --conditioning_closeness "${sigma}" \
        --num_samples "${NUM_SAMPLES}" \
        --batch_size "${BATCH_SIZE}" \
        --output_fasta "${output_dir}/2_analog_only/${proto_type}_analog_only.fasta" \
        --conditioning_output_path "${output_dir}/2_analog_only/${proto_type}_conditioning.pt"

    echo "✓ Completed: Analog only (${proto_type})"
}

# Condition 3: Analog + Targeted physicochemical (Property conditioning)
run_analog_property() {
    local prototype_file=$1
    local output_dir=$2
    local proto_type=$3
    local tau=$4
    local sigma=$5

    echo "----------------------------------------"
    echo "Condition 3: Analog + Property (${proto_type}, τ=${tau}, σ=${sigma})"
    echo "----------------------------------------"

    python "${GENERATE_SCRIPT}" AdvancedConditional \
        --analog_file "${prototype_file}" \
        --checkpoint_path "${CHECKPOINT_PATH}" \
        --length "${LENGTH_RANGE}" \
        --charge "${CHARGE_RANGE}" \
        --hydrophobicity="${HYDROPHOBICITY_RANGE}" \
        --initial_timestep "${tau}" \
        --conditioning_closeness "${sigma}" \
        --num_samples "${NUM_SAMPLES}" \
        --batch_size "${BATCH_SIZE}" \
        --output_fasta "${output_dir}/3_analog_property/${proto_type}_analog_property.fasta" \
        --conditioning_output_path "${output_dir}/3_analog_property/${proto_type}_conditioning.pt" \
        --advanced_conditioning_modes "AnalogConditional,PartialConditional"

    echo "✓ Completed: Analog + Property (${proto_type})"
}

# Condition 4: Analog + Template (dynamic motif insertion)
run_analog_template() {
    local prototype_file=$1
    local output_dir=$2
    local proto_type=$3
    local tau=$4
    local sigma=$5

    echo "----------------------------------------"
    echo "Condition 4: Analog + Template (${proto_type}, τ=${tau}, σ=${sigma})"
    echo "  Motif: ${MOTIF} (inserted at random positions)"
    echo "----------------------------------------"

    python "${GENERATE_SCRIPT}" AdvancedConditional \
        --analog_file "${prototype_file}" \
        --checkpoint_path "${CHECKPOINT_PATH}" \
        --motif "${MOTIF}" \
        --motif_seed "${MOTIF_SEED}" \
        --initial_timestep "${tau}" \
        --conditioning_closeness "${sigma}" \
        --num_samples "${NUM_SAMPLES}" \
        --batch_size "${BATCH_SIZE}" \
        --output_fasta "${output_dir}/4_analog_template/${proto_type}_analog_template.fasta" \
        --conditioning_output_path "${output_dir}/4_analog_template/${proto_type}_conditioning.pt" \
        --advanced_conditioning_modes "TemplateConditional,AnalogConditional"

    echo "✓ Completed: Analog + Template (${proto_type})"
}

# Condition 5: Analog + Template + Targeted physicochemical
run_analog_template_property() {
    local prototype_file=$1
    local output_dir=$2
    local proto_type=$3
    local tau=$4
    local sigma=$5

    echo "----------------------------------------"
    echo "Condition 5: Analog + Template + Property (${proto_type}, τ=${tau}, σ=${sigma})"
    echo "  Motif: ${MOTIF} (inserted at random positions)"
    echo "----------------------------------------"

    python "${GENERATE_SCRIPT}" AdvancedConditional \
        --analog_file "${prototype_file}" \
        --checkpoint_path "${CHECKPOINT_PATH}" \
        --motif "${MOTIF}" \
        --motif_seed "${MOTIF_SEED}" \
        --length "${LENGTH_RANGE}" \
        --charge "${CHARGE_RANGE}" \
        --hydrophobicity="${HYDROPHOBICITY_RANGE}" \
        --initial_timestep "${tau}" \
        --conditioning_closeness "${sigma}" \
        --num_samples "${NUM_SAMPLES}" \
        --batch_size "${BATCH_SIZE}" \
        --output_fasta "${output_dir}/5_analog_template_property/${proto_type}_analog_template_property.fasta" \
        --conditioning_output_path "${output_dir}/5_analog_template_property/${proto_type}_conditioning.pt" \
        --advanced_conditioning_modes "TemplateConditional,AnalogConditional,PartialConditional"

    echo "✓ Completed: Analog + Template + Property (${proto_type})"
}

run_all_conditions() {
    local prototype_file=$1
    local output_dir=$2
    local proto_type=$3
    local tau=$4
    local sigma=$5

    echo ""
    echo "========================================"
    echo "Running all conditions for ${proto_type} prototypes"
    echo "Parameters: τ=${tau}, σ=${sigma}"
    echo "========================================"
    echo ""

    copy_prototypes "${prototype_file}" "${output_dir}" "${proto_type}"
    run_analog_only "${prototype_file}" "${output_dir}" "${proto_type}" "${tau}" "${sigma}"
    run_analog_property "${prototype_file}" "${output_dir}" "${proto_type}" "${tau}" "${sigma}"
    run_analog_template "${prototype_file}" "${output_dir}" "${proto_type}" "${tau}" "${sigma}"
    run_analog_template_property "${prototype_file}" "${output_dir}" "${proto_type}" "${tau}" "${sigma}"

    echo ""
    echo "========================================"
    echo "Completed all conditions for ${proto_type}"
    echo "========================================"
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

echo "================================================================"
echo "ANALOG GENERATION FLAVORS COMPARISON"
echo "================================================================"
echo ""
echo "Model: ${CHECKPOINT_PATH}"
echo ""
echo "Baselines:"
echo "  OmegAMP-U (unconditional): ${OMEGAMP_U_SOURCE}"
echo "  OmegAMP-T (targeted): ${OMEGAMP_T_SOURCE}"
echo ""
echo "Positive prototypes: ${POSITIVE_PROTOTYPES}"
echo "  τ=${POSITIVE_TAU}, σ=${POSITIVE_SIGMA}"
echo ""
echo "Negative prototypes: ${NEGATIVE_PROTOTYPES}"
echo "  τ=${NEGATIVE_TAU}, σ=${NEGATIVE_SIGMA}"
echo ""
echo "Samples per prototype: ${NUM_SAMPLES}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Property conditioning ranges:"
echo "  Length: ${LENGTH_RANGE}"
echo "  Charge: ${CHARGE_RANGE}"
echo "  Hydrophobicity: ${HYDROPHOBICITY_RANGE}"
echo ""
echo "Template motif: ${MOTIF}"
echo "  Seed: ${MOTIF_SEED}"
echo "  (Motif will be inserted at random positions in each prototype)"
echo ""
echo "Conditions to generate:"
echo "  0. Baselines (OmegAMP-U, OmegAMP-T)"
echo "  1. Prototypes (baseline)"
echo "  2. Analog only"
echo "  3. Analog + Targeted physicochemical"
echo "  4. Analog + Template (motif insertion)"
echo "  5. Analog + Template + Targeted physicochemical"
echo ""

# Check files exist
if [ ! -f "${POSITIVE_PROTOTYPES}" ]; then
    echo "ERROR: Positive prototype file not found: ${POSITIVE_PROTOTYPES}"
    exit 1
fi

if [ ! -f "${NEGATIVE_PROTOTYPES}" ]; then
    echo "ERROR: Negative prototype file not found: ${NEGATIVE_PROTOTYPES}"
    exit 1
fi

if [ ! -f "${GENERATE_SCRIPT}" ]; then
    echo "ERROR: Generate script not found: ${GENERATE_SCRIPT}"
    exit 1
fi


if [ ! -f "${OMEGAMP_U_SOURCE}" ]; then
    echo "WARNING: OmegAMP-U source not found: ${OMEGAMP_U_SOURCE}"
fi

if [ ! -f "${OMEGAMP_T_SOURCE}" ]; then
    echo "WARNING: OmegAMP-T source not found: ${OMEGAMP_T_SOURCE}"
fi

# Prompt for confirmation
read -p "Press Enter to start generation, or Ctrl+C to cancel..."

# Record start time
start_time=$(date +%s)

# Create output directories
create_directories "${OUTPUT_DIR}"

# Copy baselines (Condition 0)
copy_omegamp_baselines "${OUTPUT_DIR}"

# Run for positive prototypes
run_all_conditions \
    "${POSITIVE_PROTOTYPES}" \
    "${OUTPUT_DIR}" \
    "positive" \
    "${POSITIVE_TAU}" \
    "${POSITIVE_SIGMA}"

# Run for negative prototypes
run_all_conditions \
    "${NEGATIVE_PROTOTYPES}" \
    "${OUTPUT_DIR}" \
    "negative" \
    "${NEGATIVE_TAU}" \
    "${NEGATIVE_SIGMA}"

# Record end time
end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo ""
echo "================================================================"
echo "GENERATION COMPLETED"
echo "================================================================"
echo "Total time: ${hours}h ${minutes}m ${seconds}s"
echo "Results saved to: ${OUTPUT_DIR}"
echo ""
echo "Output structure:"
echo "  ${OUTPUT_DIR}/"
echo "    ├── 0_omegamp_baselines/"
echo "    │   ├── OmegAMP-U.fasta  (unconditional)"
echo "    │   └── OmegAMP-T.fasta  (targeted physicochemical)"
echo "    ├── 1_prototypes/"
echo "    │   ├── positive_prototypes.fasta"
echo "    │   └── negative_prototypes.fasta"
echo "    ├── 2_analog_only/"
echo "    ├── 3_analog_property/"
echo "    ├── 4_analog_template/"
echo "    │   ├── *.fasta"
echo "    │   └── *.templates.txt  (motif positions log)"
echo "    └── 5_analog_template_property/"
echo "        ├── *.fasta"
echo "        └── *.templates.txt  (motif positions log)"
echo "================================================================"

# Save configuration
cat > "${OUTPUT_DIR}/experiment_config.txt" << EOF
Analog Generation Flavors Comparison
=====================================
Date: $(date)
Model: ${CHECKPOINT_PATH}

Baselines:
  OmegAMP-U (unconditional): ${OMEGAMP_U_SOURCE}
  OmegAMP-T (targeted): ${OMEGAMP_T_SOURCE}

Prototype files:
  Positive: ${POSITIVE_PROTOTYPES}
  Negative: ${NEGATIVE_PROTOTYPES}

Parameters:
  Positive: τ=${POSITIVE_TAU}, σ=${POSITIVE_SIGMA}
  Negative: τ=${NEGATIVE_TAU}, σ=${NEGATIVE_SIGMA}

Property conditioning:
  Length: ${LENGTH_RANGE}
  Charge: ${CHARGE_RANGE}
  Hydrophobicity: ${HYDROPHOBICITY_RANGE}

Template conditioning:
  Motif: ${MOTIF}
  Seed: ${MOTIF_SEED}
  Method: Dynamic insertion at random positions per prototype
  Note: Each prototype receives a length-matched template with the
        motif inserted at a reproducible random position.

Sampling:
  Samples per prototype: ${NUM_SAMPLES}
  Batch size: ${BATCH_SIZE}

Conditions:
  0. Baselines
     - OmegAMP-U: Unconditional generation
     - OmegAMP-T: Targeted physicochemical generation
  1. Prototypes (baseline)
  2. Analog only (AnalogConditional)
  3. Analog + Targeted physicochemical (AnalogConditional + PartialConditional)
  4. Analog + Template (AnalogConditional + TemplateConditional)
  5. Analog + Template + Targeted (AnalogConditional + TemplateConditional + PartialConditional)

Duration: ${hours}h ${minutes}m ${seconds}s
EOF

echo ""
echo "Configuration saved to: ${OUTPUT_DIR}/experiment_config.txt"
echo ""