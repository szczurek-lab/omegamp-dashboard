#!/usr/bin/env bash
# Inactive-to-active wet-lab analog generation. Six inactive prototypes, three
# strategies per prototype:
#
#   AU (broad sampling):  analog unconditional at tau=0.25, sigma=0.25
#   AT (closeness sweep + property-targeted):
#                         analog unconditional, sigma in [0, 1] at tau=0.1
#                         analog targeted, standard property ranges, tau=0.25
#   AP (prototype-derived): analog prototype-derived with the parent as
#                           prototype source, tau=0.25, sigma=0.5
#
# 50,000 candidates per cell. Run under the OmegAMP env.
# Usage: bash scripts/reproduce/wetlab_02_generate_analog.sh [--mini] [--only NAME]
set -euo pipefail

source config.sh

declare -A PROTOTYPES=(
    [Mammutin-1]="WMTIHALKLSLSFKL"
    [GQ20]="GQLNKFIKKAQRKFHEKFAK"
    [BoCo1]="NKIKFINKYVKKVQLKKILVKS"
    [DeNo1047]="ALPSIIKGLLKKL"
    [As-CATH4-6L]="IGKLFKRIVERIKRFLRVLLRILR"
    [OP-145-TII4]="LGLFKLLLRLILKGFKL"
)

LENGTH="5:30"
CHARGE="2:10"
HYDRO="-0.5:0.8"

CLOSENESS_VALUES=(0.0 0.25 0.5 0.75 1.0)
SAMPLES=50000
ONLY=""

while [ $# -gt 0 ]; do
    case "$1" in
        --mini) SAMPLES=20; CLOSENESS_VALUES=(0.0 1.0) ;;
        --only) ONLY="$2"; shift ;;
    esac
    shift
done

OUT="${DATA_DIR}/wetlab/analog"

for name in Mammutin-1 GQ20 BoCo1 DeNo1047 As-CATH4-6L OP-145-TII4; do
    [ -n "${ONLY}" ] && [ "${ONLY}" != "${name}" ] && continue
    parent="${PROTOTYPES[$name]}"
    proto_dir="${OUT}/${name}"
    mkdir -p "${proto_dir}"

    parent_fasta="${proto_dir}/parent.fasta"
    parent_rep="${proto_dir}/parent_repeated.fasta"
    printf '>%s\n%s\n' "${name}" "${parent}" > "${parent_fasta}"
    [ -s "${parent_rep}" ] || repeat_fasta "${parent_fasta}" "${SAMPLES}" "${parent_rep}"

    # AT: analog unconditional, sigma sweep at tau=0.1
    for sigma in "${CLOSENESS_VALUES[@]}"; do
        cell="${proto_dir}/unconditional/sigma${sigma}"
        [ -s "${cell}/samples.fasta" ] && continue
        mkdir -p "${cell}"
        cd "${OMEGAMP_DIR}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog unconditional \
            --analog_sequences "${parent_rep}" \
            --tau 0.1 --sigma "${sigma}" \
            --checkpoint_path "${CHECKPOINT}" \
            --num_samples "${SAMPLES}" --batch_size 1000 \
            --output_fasta "${cell}/samples.fasta" \
            --conditioning_output_path "${cell}/conditioning.pt"
    done

    # AT: analog targeted, standard property ranges
    cell="${proto_dir}/targeted"
    if [ ! -s "${cell}/samples.fasta" ]; then
        mkdir -p "${cell}"
        cd "${OMEGAMP_DIR}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog targeted \
            --analog_sequences "${parent_rep}" \
            --tau 0.25 --sigma 0.25 \
            --length="${LENGTH}" --charge="${CHARGE}" --hydrophobicity="${HYDRO}" \
            --checkpoint_path "${CHECKPOINT}" \
            --num_samples "${SAMPLES}" --batch_size 1000 \
            --output_fasta "${cell}/samples.fasta" \
            --conditioning_output_path "${cell}/conditioning.pt"
    fi

    # AU: analog unconditional. Broad sampling from the AMP space, no property
    # or template constraints; tau and sigma match the original "analog-random"
    # campaign's timing.
    cell="${proto_dir}/unconditional-au"
    if [ ! -s "${cell}/samples.fasta" ]; then
        mkdir -p "${cell}"
        cd "${OMEGAMP_DIR}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog unconditional \
            --analog_sequences "${parent_rep}" \
            --tau 0.25 --sigma 0.25 \
            --checkpoint_path "${CHECKPOINT}" \
            --num_samples "${SAMPLES}" --batch_size 1000 \
            --output_fasta "${cell}/samples.fasta" \
            --conditioning_output_path "${cell}/conditioning.pt"
    fi

    # AP: analog prototype-derived. Parent peptide provides both the analog
    # seed and the prototype conditioning.
    cell="${proto_dir}/prototype-derived"
    if [ ! -s "${cell}/samples.fasta" ]; then
        mkdir -p "${cell}"
        cd "${OMEGAMP_DIR}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog prototype-derived \
            --analog_sequences "${parent_rep}" \
            --prototype_sequences "${parent_rep}" \
            --tau 0.25 --sigma 0.5 \
            --checkpoint_path "${CHECKPOINT}" \
            --num_samples "${SAMPLES}" --batch_size 1000 \
            --output_fasta "${cell}/samples.fasta" \
            --conditioning_output_path "${cell}/conditioning.pt"
    fi
done

echo "Done. $(find "${OUT}" -name 'samples.fasta' | wc -l) FASTA files in ${OUT}"