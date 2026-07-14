#!/usr/bin/env bash
# Compute mean MMseqs2 similarity-to-prototype for every figure 3 flavors cell.
# De novo flavors are compared against the positive (active AMPs) set as a
# reference; analog flavors against their own prototype set.
# Needs the amp-modelling env (has mmseqs2).
# Usage: bash scripts/reproduce/compute_similarity_figure3_flavors.sh [--only positive|negative]
set -euo pipefail

source config.sh

ONLY=both
while [ $# -gt 0 ]; do
    case "$1" in
        --only) ONLY="$2"; shift ;;
    esac
    shift
done

SCRIPT="${DASHBOARD_DIR}/scripts/mean_mmseqs_score.py"
FLAVORS_DIR="${DATA_DIR}/figure3/flavors"

# Each cell's similarity is computed in its own tempdir so mmseqs scratch
# files don't collide between cells.
score_one() {
    local fasta="$1" database="$2" out_csv="$3"
    [ -s "${out_csv}" ] && return 0
    mkdir -p "$(dirname "${out_csv}")"
    local tmp; tmp=$(mktemp -d)
    ( cd "${tmp}" && \
      python "${SCRIPT}" --query "${fasta}" --database "${database}" --output "${out_csv}" )
    rm -rf "${tmp}"
}

# De novo flavors (no set axis): compare vs positive prototypes by default.
if [ "${ONLY}" != "negative" ]; then
    for flavor_stem in \
        "omegamp_du/OmegAMP-U" \
        "omegamp_dt/OmegAMP-T" \
        "omegamp_dp/denovo_subset"; do
        flavor_dir="${flavor_stem%/*}"
        stem="${flavor_stem##*/}"
        fasta="${FLAVORS_DIR}/${flavor_dir}/${stem}.fasta"
        out_csv="${FLAVORS_DIR}/${flavor_dir}/similarity/${stem}_vs_positive_prototypes.csv"
        [ -s "${fasta}" ] && score_one "${fasta}" "${POSITIVE_PROTOTYPES}" "${out_csv}"
    done
fi

# Analog flavors: each set compared against its own prototype set.
for set_name in positive negative; do
    [ "${ONLY}" != "both" ] && [ "${ONLY}" != "${set_name}" ] && continue
    if [ "${set_name}" = "positive" ]; then
        database="${POSITIVE_PROTOTYPES}"
    else
        database="${NEGATIVE_PROTOTYPES}"
    fi
    for flavor_stem in \
        "omegamp_au/${set_name}_analog_only" \
        "omegamp_at/${set_name}_analog_property" \
        "omegamp_am/${set_name}_analog_template" \
        "omegamp_amt/${set_name}_analog_template_property" \
        "omegamp_ap/${set_name}_analog_subset"; do
        flavor_dir="${flavor_stem%/*}"
        stem="${flavor_stem##*/}"
        fasta="${FLAVORS_DIR}/${flavor_dir}/${stem}.fasta"
        out_csv="${FLAVORS_DIR}/${flavor_dir}/similarity/${stem}_vs_${set_name}_prototypes.csv"
        [ -s "${fasta}" ] && score_one "${fasta}" "${database}" "${out_csv}"
    done
done

echo "Done. $(find "${FLAVORS_DIR}" -name '*_vs_*_prototypes.csv' | wc -l) similarity files."
