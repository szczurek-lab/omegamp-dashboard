#!/usr/bin/env bash
# OmegAMP de-novo benchmark sets for figure 2 Panel H (and figure 2/S1 novelty).
# Three generation modes, written to data/figure2/benchmarks/:
#   unconditional.fasta          OmegAMP-DU  de-novo unconditional
#   target-physicochemical.fasta OmegAMP-DT  de-novo targeted (length/charge/hydrophobicity)
#   population-guided.fasta      OmegAMP-DP  de-novo prototype-derived (positive prototypes)
#
# These are the same generation modes as the figure 3 de-novo flavors; here they
# are emitted as standalone benchmark FASTAs. They are used for sequence-property
# and novelty panels only (no APEX scoring needed). The competitor FASTAs
# (amp-gan, amp-diffusion, ...), species MIC tables, and the training set are
# external inputs obtained from the data deposit -- see README.
#
# Run under the OmegAMP env (omegamp-inference).
# Usage: bash scripts/reproduce/gen_figure2_benchmarks.sh [--mini]
set -euo pipefail

source config.sh

if [ "${1:-}" = "--mini" ]; then
    NUM_SAMPLES=100
else
    # Property distributions are insensitive to N; the deposited benchmark sets
    # used larger draws (~50k-150k). Override with NUM_SAMPLES=... if matching exactly.
    NUM_SAMPLES="${NUM_SAMPLES:-50000}"
fi

LENGTH="5:30"
CHARGE="2:10"
HYDRO="-0.5:0.8"

OUT="${DATA_DIR}/figure2/benchmarks"
mkdir -p "${OUT}"

cd "${OMEGAMP_DIR}"

# OmegAMP-DU -- de-novo unconditional
fasta="${OUT}/unconditional.fasta"
if [ ! -s "${fasta}" ]; then
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo unconditional \
        --checkpoint_path "${CHECKPOINT}" \
        --num_samples "${NUM_SAMPLES}" --batch_size 256 \
        --output_fasta "${fasta}" \
        --conditioning_output_path "${OUT}/unconditional_conditioning.pt"
fi

# OmegAMP-DT -- de-novo targeted (physicochemical)
fasta="${OUT}/target-physicochemical.fasta"
if [ ! -s "${fasta}" ]; then
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
        --checkpoint_path "${CHECKPOINT}" \
        --length="${LENGTH}" --charge="${CHARGE}" --hydrophobicity="${HYDRO}" \
        --num_samples "${NUM_SAMPLES}" --batch_size 256 \
        --output_fasta "${fasta}" \
        --conditioning_output_path "${OUT}/target-physicochemical_conditioning.pt"
fi

# OmegAMP-DP -- de-novo prototype-derived (positive prototypes)
fasta="${OUT}/population-guided.fasta"
if [ ! -s "${fasta}" ]; then
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo prototype-derived \
        --checkpoint_path "${CHECKPOINT}" \
        --prototype_sequences "${POSITIVE_PROTOTYPES}" \
        --sigma "${OPTIMAL_POSITIVE_SIGMA}" \
        --num_samples "${NUM_SAMPLES}" --batch_size 256 \
        --output_fasta "${fasta}" \
        --conditioning_output_path "${OUT}/population-guided_conditioning.pt"
fi

echo "Done. Benchmark FASTAs in ${OUT}:"
for f in unconditional target-physicochemical population-guided; do
    printf "  %-25s %s seqs\n" "${f}.fasta" "$(grep -c '^>' "${OUT}/${f}.fasta" 2>/dev/null || echo 0)"
done
