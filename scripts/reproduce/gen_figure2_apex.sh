#!/usr/bin/env bash
# =============================================================================
# Run APEX MIC predictions for Figure 2 grid-sweep samples.
#
# Requires gen_figure2_conditioning.sh to have been run first.
#
# Reads:  data/figure2/conditioning/grid-sweep-charge-hydrophob/c{C}_h{H}/samples.fasta
# Writes: data/figure2/apex/grid-sweep-charge-hydrophob/c{C}_h{H}-apex-predictions-min.tsv
#
# The "-min" suffix is a naming convention; the TSV contains all-strain
# predictions from APEX and the notebook takes min MIC per Sequence_id itself.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/config.sh"
check_config
check_apex

COND_DIR="${DATA_DIR}/figure2/conditioning/grid-sweep-charge-hydrophob"
APEX_OUT="${DATA_DIR}/figure2/apex/grid-sweep-charge-hydrophob"

CHARGE_GRID=(0 2 4 6 8 10 12)
HYDROPHOB_GRID=("-0.75" "-0.50" "-0.25" "0.0" "0.25" "0.5" "0.75")

total=$(( ${#CHARGE_GRID[@]} * ${#HYDROPHOB_GRID[@]} ))
done_count=0
skipped=0

echo "================================================================"
echo "Figure 2 — APEX predictions (Charge × Hydrophobicity grid)"
echo "APEX:   ${APEX_DIR}"
echo "Input:  ${COND_DIR}"
echo "Output: ${APEX_OUT}"
echo "Grid points: ${total}"
echo "================================================================"

mkdir -p "${APEX_OUT}"

for C in "${CHARGE_GRID[@]}"; do
    for H in "${HYDROPHOB_GRID[@]}"; do
        HC=$(fmt_hydro "${H}")
        DN="c${C}_h${HC}"
        INPUT="${COND_DIR}/${DN}/samples.fasta"
        OUTPUT="${APEX_OUT}/${DN}-apex-predictions-min.tsv"

        if [ ! -f "${INPUT}" ]; then
            echo "  [${DN}] SKIP — samples.fasta not found (run gen_figure2_conditioning.sh first)"
            (( skipped++ )) || true
            continue
        fi

        if [ -f "${OUTPUT}" ] && [ -s "${OUTPUT}" ]; then
            echo "  [${DN}] skip (predictions already exist)"
            (( done_count++ )) || true
            continue
        fi

        echo "  [${DN}] C=${C}, H=${H}"
        run_apex "${INPUT}" "${OUTPUT}"
        (( done_count++ )) || true
    done
done

echo ""
echo "================================================================"
echo "APEX predictions complete."
echo "  Generated / already present: ${done_count}"
echo "  Skipped (missing input):     ${skipped}"
echo "  Output: ${APEX_OUT}"
echo "================================================================"
