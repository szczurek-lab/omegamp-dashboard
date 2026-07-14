#!/usr/bin/env bash
# De-novo property-conditioned samples for figure 2.
# Usage: bash scripts/reproduce/gen_figure2_conditioning.sh [--mini]
set -euo pipefail

source config.sh

if [ "${1:-}" = "--mini" ]; then
    NUM_SAMPLES=10
else
    NUM_SAMPLES=500
fi

OUT="${DATA_DIR}/figure2/conditioning"

COMMON_ARGS=(
    --checkpoint_path "${CHECKPOINT}"
    --num_samples "${NUM_SAMPLES}"
    --batch_size "${NUM_SAMPLES}"
)

# Hydra resolves its config relative to the OmegAMP repo.
cd "${OMEGAMP_DIR}"

# ---------------------------------------------------------------------------
# Single-property sweeps
# ---------------------------------------------------------------------------
for c in 0 2 4 6 8 10; do
    dir="${OUT}/single-property/charge-${c}"
    mkdir -p "${dir}"
    [ -s "${dir}/samples.fasta" ] && continue
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
        --charge "${c}" \
        --output_fasta "${dir}/samples.fasta" \
        --conditioning_output_path "${dir}/conditioning.pt" \
        "${COMMON_ARGS[@]}"
done

for l in 10 15 20 25 30; do
    dir="${OUT}/single-property/length-${l}"
    mkdir -p "${dir}"
    [ -s "${dir}/samples.fasta" ] && continue
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
        --length "${l}" \
        --output_fasta "${dir}/samples.fasta" \
        --conditioning_output_path "${dir}/conditioning.pt" \
        "${COMMON_ARGS[@]}"
done

for h in -0.5 -0.2 0.0 0.2 0.4 0.6 0.8; do
    h_tag=$(echo "${h}" | tr '.-' 'pm')
    dir="${OUT}/single-property/hydrophob-${h_tag}"
    mkdir -p "${dir}"
    [ -s "${dir}/samples.fasta" ] && continue
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
        --hydrophobicity "${h}" \
        --output_fasta "${dir}/samples.fasta" \
        --conditioning_output_path "${dir}/conditioning.pt" \
        "${COMMON_ARGS[@]}"
done

# ---------------------------------------------------------------------------
# Multi-property conditions
# ---------------------------------------------------------------------------
MULTI=(
    low_charge_short:2:10:-0.2
    moderate_low:4:15:0.0
    moderate:6:20:0.3
    high_moderate:8:25:0.4
    high_long:10:30:0.6
)

for entry in "${MULTI[@]}"; do
    IFS=: read -r name c l h <<< "${entry}"
    dir="${OUT}/multi-property/${name}"
    mkdir -p "${dir}"
    [ -s "${dir}/samples.fasta" ] && continue
    CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
        --charge "${c}" \
        --length "${l}" \
        --hydrophobicity "${h}" \
        --output_fasta "${dir}/samples.fasta" \
        --conditioning_output_path "${dir}/conditioning.pt" \
        "${COMMON_ARGS[@]}"
done

# ---------------------------------------------------------------------------
# 2D grid sweeps
# ---------------------------------------------------------------------------
for c in 0 2 4 6 8 10 12; do
    for h in -0.75 -0.50 -0.25 0.0 0.25 0.5 0.75; do
        h_tag=$(echo "${h}" | tr '.-' 'pm')
        dir="${OUT}/grid-sweep-charge-hydrophob/c${c}_h${h_tag}"
        mkdir -p "${dir}"
        [ -s "${dir}/samples.fasta" ] && continue
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
            --charge "${c}" \
            --hydrophobicity "${h}" \
            --output_fasta "${dir}/samples.fasta" \
            --conditioning_output_path "${dir}/conditioning.pt" \
            "${COMMON_ARGS[@]}"
    done
done

for c in 0 2 4 6 8 10 12; do
    for l in 5 10 15 20 25 30 35; do
        dir="${OUT}/grid-sweep-charge-length/c${c}_l${l}"
        mkdir -p "${dir}"
        [ -s "${dir}/samples.fasta" ] && continue
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
            --charge "${c}" \
            --length "${l}" \
            --output_fasta "${dir}/samples.fasta" \
            --conditioning_output_path "${dir}/conditioning.pt" \
            "${COMMON_ARGS[@]}"
    done
done

for h in -0.75 -0.50 -0.25 0.0 0.25 0.5 0.75; do
    h_tag=$(echo "${h}" | tr '.-' 'pm')
    for l in 5 10 15 20 25 30 35; do
        dir="${OUT}/grid-sweep-hydrophob-length/h${h_tag}_l${l}"
        mkdir -p "${dir}"
        [ -s "${dir}/samples.fasta" ] && continue
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" de-novo targeted \
            --hydrophobicity "${h}" \
            --length "${l}" \
            --output_fasta "${dir}/samples.fasta" \
            --conditioning_output_path "${dir}/conditioning.pt" \
            "${COMMON_ARGS[@]}"
    done
done

echo "Done. $(find "${OUT}" -name samples.fasta | wc -l) FASTA files in ${OUT}"