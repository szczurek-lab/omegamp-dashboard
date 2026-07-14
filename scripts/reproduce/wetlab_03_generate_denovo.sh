#!/usr/bin/env bash
# De novo wet-lab generation:
#   OmegAMP-P: de-novo prototype-derived, conditioned on curated AMPs
#   OmegAMP-T: de-novo targeted, standard property ranges
# 50,000 candidates each. Run under the OmegAMP env.
# Usage: bash scripts/reproduce/wetlab_03_generate_denovo.sh [--mini]
set -euo pipefail

source config.sh

LENGTH="10:30"
CHARGE="2:10"
HYDRO="-0.5:0.8"
SAMPLES=50000

while [ $# -gt 0 ]; do
    case "$1" in
        --mini) SAMPLES=100 ;;
    esac
    shift
done

OUT="${DATA_DIR}/wetlab/denovo"
mkdir -p "${OUT}"

# OmegAMP-P: de novo prototype-derived. The new API uses the prototype-FASTA
# entry count as the output count, so repeat the curated AMPs to reach SAMPLES.
cell="${OUT}/prototype-derived"
if [ ! -s "${cell}/samples.fasta" ]; then
    mkdir -p "${cell}"
    amps_rep="${cell}/amps_repeated.fasta"
    n_proto=$(grep -c '^>' "${POSITIVE_PROTOTYPES}")
    repeats=$(( (SAMPLES + n_proto - 1) / n_proto ))
    repeat_fasta "${POSITIVE_PROTOTYPES}" "${repeats}" "${amps_rep}"
    cd "${OMEGAMP_DIR}"
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo prototype-derived \
        --prototype_sequences "${amps_rep}" \
        --sigma 0.0 \
        --checkpoint_path "${CHECKPOINT}" \
        --num_samples "$(grep -c '^>' "${amps_rep}")" --batch_size 256 \
        --output_fasta "${cell}/samples.fasta" \
        --conditioning_output_path "${cell}/conditioning.pt"
fi

# OmegAMP-T: de novo targeted, standard property ranges
cell="${OUT}/targeted"
if [ ! -s "${cell}/samples.fasta" ]; then
    mkdir -p "${cell}"
    cd "${OMEGAMP_DIR}"
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
        --length="${LENGTH}" --charge="${CHARGE}" --hydrophobicity="${HYDRO}" \
        --checkpoint_path "${CHECKPOINT}" \
        --num_samples "${SAMPLES}" --batch_size 256 \
        --output_fasta "${cell}/samples.fasta" \
        --conditioning_output_path "${cell}/conditioning.pt"
fi

echo "Done. $(find "${OUT}" -name 'samples.fasta' | wc -l) FASTA files in ${OUT}"