#!/usr/bin/env bash
# Compute mean MMseqs2 similarity-to-prototype for every figure 3 sweep cell.
# Needs the amp-modelling env (has mmseqs2).
# Usage: bash scripts/reproduce/compute_similarity_figure3_sweep.sh [--only positive|negative]
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

# Each cell's similarity is computed in its own tempdir so mmseqs scratch
# files (result_*.tsv, tmp/) don't collide between cells.
score_one() {
    local fasta="$1" database="$2" out_csv="$3"
    [ -s "${out_csv}" ] && return 0
    mkdir -p "$(dirname "${out_csv}")"
    local tmp; tmp=$(mktemp -d)
    ( cd "${tmp}" && \
      python "${SCRIPT}" --query "${fasta}" --database "${database}" --output "${out_csv}" )
    rm -rf "${tmp}"
}

for set_name in positive negative; do
    [ "${ONLY}" != "both" ] && [ "${ONLY}" != "${set_name}" ] && continue
    if [ "${set_name}" = "positive" ]; then
        database="${POSITIVE_PROTOTYPES}"
    else
        database="${NEGATIVE_PROTOTYPES}"
    fi

    seq_dir="${DATA_DIR}/figure3/sweep/${set_name}/sequences"
    sim_dir="${DATA_DIR}/figure3/sweep/${set_name}/similarity"

    for fasta in "${seq_dir}"/tau*.fasta; do
        cell=$(basename "${fasta}" .fasta)
        out_csv="${sim_dir}/${cell}_vs_${set_name}_prototypes.csv"
        score_one "${fasta}" "${database}" "${out_csv}"
    done
done

echo "Done. $(find "${DATA_DIR}/figure3/sweep" -name '*_vs_*_prototypes.csv' | wc -l) similarity files."
