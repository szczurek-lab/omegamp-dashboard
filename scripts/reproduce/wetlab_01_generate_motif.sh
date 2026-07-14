#!/usr/bin/env bash
# Motif-guided wet-lab generation.
# LPS-binding motifs: analog-motif prototype-derived, sweep sigma at tau=0.25.
# bZIP: motif targeted with property ranges, single configuration.
# Filtering and selection happen in later steps. Run under the OmegAMP env.
# Usage: bash scripts/reproduce/wetlab_01_generate_motif.sh [--mini] [--only MOTIF]
set -euo pipefail

source config.sh

# LPS-binding motif templates and parent peptides.
declare -A LPS_TEMPLATES=(
    [cecropin]="--------------KK-----I---I-----"
    [sarcotoxin]="-W-KK----------------------------------"
    [pa4]="-------KII--P--K-LL-A------------"
    [LG21]="-------------GWKRKRFG"
)
declare -A LPS_PARENTS=(
    [cecropin]="SWLSKTAKKLENSAKKRISEGIAIAIQGGPR"
    [sarcotoxin]="GWLKKIGKKIERVGQHTRDATIQGLGIAQQAANVAATAR"
    [pa4]="GFFALIPKIISSPLFKTLLSAVGSALSSSGGQE"
    [LG21]="LLPIVGNLLKSLLGWKRKRFG"
)

BZIP_TEMPLATE="---L---V--L---N--L---V--L---V-"

LPS_TAU=0.25
LPS_SIGMA_VALUES=(0.0 0.25 0.5 0.75 1.0)
SAMPLES=50000
ONLY=""

while [ $# -gt 0 ]; do
    case "$1" in
        --mini)
            SAMPLES=20
            LPS_SIGMA_VALUES=(0.0 1.0)
            ;;
        --only) ONLY="$2"; shift ;;
    esac
    shift
done

OUT="${DATA_DIR}/wetlab/motif"

# ---------------------------------------------------------------------------
# LPS-binding motifs: analog-motif prototype-derived, sweep sigma.
# ---------------------------------------------------------------------------
for motif in cecropin sarcotoxin pa4 LG21; do
    [ -n "${ONLY}" ] && [ "${ONLY}" != "${motif}" ] && continue
    template="${LPS_TEMPLATES[$motif]}"
    parent="${LPS_PARENTS[$motif]}"

    for sigma in "${LPS_SIGMA_VALUES[@]}"; do
        cell="${OUT}/${motif}/sigma${sigma}"
        [ -s "${cell}/samples.fasta" ] && continue
        mkdir -p "${cell}"

        printf '>%s\n%s\n' "${motif}" "${template}" > "${cell}/motif.fasta"
        printf '>%s\n%s\n' "${motif}" "${parent}"   > "${cell}/parent.fasta"
        repeat_fasta "${cell}/motif.fasta"  "${SAMPLES}" "${cell}/motif_repeated.fasta"
        repeat_fasta "${cell}/parent.fasta" "${SAMPLES}" "${cell}/parent_repeated.fasta"

        cd "${OMEGAMP_DIR}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" analog-motif prototype-derived \
            --motif_sequences "${cell}/motif_repeated.fasta" \
            --analog_sequences "${cell}/parent_repeated.fasta" \
            --prototype_sequences "${cell}/parent_repeated.fasta" \
            --tau "${LPS_TAU}" --sigma "${sigma}" --guidance_strength 1.0 \
            --checkpoint_path "${CHECKPOINT}" \
            --num_samples "${SAMPLES}" --batch_size 1000 \
            --output_fasta "${cell}/samples.fasta" \
            --conditioning_output_path "${cell}/conditioning.pt"
    done
done

# ---------------------------------------------------------------------------
# bZIP: motif targeted with property ranges, single configuration.
# ---------------------------------------------------------------------------
if [ -z "${ONLY}" ] || [ "${ONLY}" = "bZIP" ]; then
    cell="${OUT}/bZIP"
    if [ ! -s "${cell}/samples.fasta" ]; then
        mkdir -p "${cell}"
        printf '>bZIP\n%s\n' "${BZIP_TEMPLATE}" > "${cell}/motif.fasta"
        repeat_fasta "${cell}/motif.fasta" "${SAMPLES}" "${cell}/motif_repeated.fasta"

        cd "${OMEGAMP_DIR}"
        CUDA_VISIBLE_DEVICES=${GPU} python "${GENERATE_SCRIPT}" motif targeted \
            --motif_sequences "${cell}/motif_repeated.fasta" \
            --charge=2:10 --hydrophobicity=-0.5:0.8 --length="${#BZIP_TEMPLATE}" \
            --guidance_strength 10.0 \
            --checkpoint_path "${CHECKPOINT}" \
            --num_samples "${SAMPLES}" --batch_size 1000 \
            --output_fasta "${cell}/samples.fasta" \
            --conditioning_output_path "${cell}/conditioning.pt"
    fi
fi

echo "Done. $(find "${OUT}" -name 'samples.fasta' | wc -l) FASTA files in ${OUT}"