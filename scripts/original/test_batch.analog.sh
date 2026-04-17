#!/bin/bash

# Quick test script for batch analog generation functionality
# Use this to verify everything works before running the full sweep

set -e

echo "=========================================="
echo "Testing Batch Analog Generation"
echo "=========================================="

# Configuration
CHECKPOINT_PATH="models/generative_model.ckpt"
TEST_ANALOGS="data/activity-data/curated-AMPs_samples.fasta"
OUTPUT_DIR="results/batch-analog-test"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Test parameters
TAU=250
SIGMA=0.5
NUM_SAMPLES=100  # PER PROTOTYPE (not total)
BATCH_SIZE=50

echo ""
echo "Test configuration:"
echo "  Checkpoint: ${CHECKPOINT_PATH}"
echo "  Prototypes: ${TEST_ANALOGS}"
echo "  τ: ${TAU}"
echo "  σ: ${SIGMA}"
echo "  Samples per prototype: ${NUM_SAMPLES}"
echo "  Batch size: ${BATCH_SIZE}"
echo ""

# Check if test file exists
if [ ! -f "${TEST_ANALOGS}" ]; then
    echo "ERROR: Test prototypes file not found: ${TEST_ANALOGS}"
    echo ""
    echo "Please run first:"
    echo "  python scripts/analysis/create_prototype_samples.py --n_samples 10"
    exit 1
fi

# Count prototypes
N_PROTOTYPES=$(grep -c "^>" "${TEST_ANALOGS}" || echo "0")
TOTAL_SAMPLES=$((N_PROTOTYPES * NUM_SAMPLES))

echo "Found ${N_PROTOTYPES} prototypes in test file"
echo "Will generate ${TOTAL_SAMPLES} total sequences (${NUM_SAMPLES} per prototype)"
echo ""

# Run test
echo "Running batch analog generation..."
python project/scripts/inference/generate_samples.py AnalogConditional \
    --checkpoint_path "${CHECKPOINT_PATH}" \
    --analog_file "${TEST_ANALOGS}" \
    --initial_timestep "${TAU}" \
    --conditioning_closeness "${SIGMA}" \
    --num_samples "${NUM_SAMPLES}" \
    --batch_size "${BATCH_SIZE}" \
    --output_fasta "${OUTPUT_DIR}/test_sequences.fasta" \
    --conditioning_output_path "${OUTPUT_DIR}/test_conditioning.pt" \
    --metadata_output_path "${OUTPUT_DIR}/test_metadata.pt"

echo ""
echo "✓ Test completed successfully!"
echo ""

# Display results
echo "Checking outputs..."
if [ -f "${OUTPUT_DIR}/test_sequences.fasta" ]; then
    num_seqs=$(grep -c "^>" "${OUTPUT_DIR}/test_sequences.fasta" || true)
    echo "  ✓ Generated ${num_seqs} sequences"

    if [ "${num_seqs}" -eq "${TOTAL_SAMPLES}" ]; then
        echo "    ✓ Count matches expected: ${TOTAL_SAMPLES}"
    else
        echo "    ⚠ WARNING: Expected ${TOTAL_SAMPLES} but got ${num_seqs}"
    fi
fi

if [ -f "${OUTPUT_DIR}/test_metadata.pt" ]; then
    echo "  ✓ Metadata saved"
    python -c "
import torch
meta = torch.load('${OUTPUT_DIR}/test_metadata.pt', weights_only=False)
print(f'    - Number of prototypes: {meta[\"n_analogs\"]}')
print(f'    - Samples per prototype: {meta[\"samples_per_analog\"]}')
print(f'    - Total samples: {len(meta[\"analog_indices\"])}')
if meta['n_analogs'] != ${N_PROTOTYPES}:
    print(f'    ⚠ WARNING: Expected {${N_PROTOTYPES}} prototypes')
if meta['samples_per_analog'] != ${NUM_SAMPLES}:
    print(f'    ⚠ WARNING: Expected {${NUM_SAMPLES}} samples per prototype')
"
fi

if [ -f "${OUTPUT_DIR}/test_conditioning.pt" ]; then
    echo "  ✓ Conditioning saved"
    python -c "
import torch
cond = torch.load('${OUTPUT_DIR}/test_conditioning.pt', weights_only=False)
print(f'    - Conditioning shape: {cond[0].shape}')
"
fi

echo ""
echo "=========================================="
echo "Test passed! Ready for parameter sweep."
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Review test outputs in ${OUTPUT_DIR}/"
echo "  2. Run full sweep: bash scripts/analysis/run_parameter_sweep.sh"