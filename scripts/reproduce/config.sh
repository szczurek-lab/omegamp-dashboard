#!/usr/bin/env bash
# =============================================================================
# Shared configuration for OmegAMP dataset generation scripts.
#
# Required environment variables (no defaults — must be set before running):
#
#   OMEGAMP_DIR   Path to the generative-modelling-amp repository
#                 e.g. export OMEGAMP_DIR=/home/user/code/generative-modelling-amp
#
#   APEX_DIR      Path to the APEX MIC prediction tool
#                 e.g. export APEX_DIR=/raid/battleamp_root/.../apex
#
# Optional overrides:
#
#   CHECKPOINT    Model checkpoint (defaults to ${OMEGAMP_DIR}/data/generative_model.ckpt)
#   GPU           CUDA device index (defaults to 0)
#
# Example usage:
#   export OMEGAMP_DIR=/home/user/code/generative-modelling-amp
#   export APEX_DIR=/raid/battleamp_root/battleamp-benchmark/battleamp/models/apex
#   bash scripts/reproduce/gen_figure2_conditioning.sh
# =============================================================================

# Dashboard repository root — auto-detected two levels above this file
_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DASHBOARD_DIR="${DASHBOARD_DIR:-$( cd "${_SCRIPT_DIR}/../.." && pwd )}"

# External paths — must be set by the user
OMEGAMP_DIR="${OMEGAMP_DIR:-}"
APEX_DIR="${APEX_DIR:-}"

# Model checkpoint
CHECKPOINT="${CHECKPOINT:-${OMEGAMP_DIR}/data/generative_model.ckpt}"

# All generated data lands here (mirrors what the notebooks read)
DATA_DIR="${DASHBOARD_DIR}/data"

# Generation entry-point scripts
GENERATE_SCRIPT="${OMEGAMP_DIR}/project/scripts/inference/generate_samples.py"
GENERATE_ANALOG_SCRIPT="${OMEGAMP_DIR}/project/scripts/inference/generate_analog_samples.py"

# Prototype FASTA files used as analog seeds
POSITIVE_PROTOTYPES="${OMEGAMP_DIR}/data/activity-data/curated-AMPs_samples.fasta"
NEGATIVE_PROTOTYPES="${OMEGAMP_DIR}/data/activity-data/curated-Non-AMPs_samples.fasta"

# GPU device index
GPU="${GPU:-0}"

# =============================================================================
# Sanity-check function — call at the top of each script
# =============================================================================
check_config() {
    local ok=1

    if [ -z "${OMEGAMP_DIR}" ]; then
        echo "ERROR: OMEGAMP_DIR is not set." >&2
        echo "  export OMEGAMP_DIR=/path/to/generative-modelling-amp" >&2
        ok=0
    elif [ ! -d "${OMEGAMP_DIR}" ]; then
        echo "ERROR: OMEGAMP_DIR not found: ${OMEGAMP_DIR}" >&2; ok=0
    fi

    if [ ! -f "${CHECKPOINT}" ]; then
        echo "ERROR: model checkpoint not found: ${CHECKPOINT}" >&2
        echo "  Override: export CHECKPOINT=/path/to/generative_model.ckpt" >&2
        ok=0
    fi
    if [ ! -f "${GENERATE_SCRIPT}" ]; then
        echo "ERROR: generate_samples.py not found: ${GENERATE_SCRIPT}" >&2; ok=0
    fi
    if [ ! -f "${GENERATE_ANALOG_SCRIPT}" ]; then
        echo "ERROR: generate_analog_samples.py not found: ${GENERATE_ANALOG_SCRIPT}" >&2; ok=0
    fi

    [ $ok -eq 1 ] || { echo "Aborting — fix the paths above." >&2; exit 1; }
}

check_apex() {
    if [ -z "${APEX_DIR}" ]; then
        echo "ERROR: APEX_DIR is not set." >&2
        echo "  export APEX_DIR=/path/to/apex" >&2
        exit 1
    fi
    if [ ! -d "${APEX_DIR}" ]; then
        echo "ERROR: APEX_DIR not found: ${APEX_DIR}" >&2; exit 1
    fi
    if [ ! -f "${APEX_DIR}/predict.py" ]; then
        echo "ERROR: predict.py not found in ${APEX_DIR}" >&2; exit 1
    fi
}

# Run APEX predictions for a single FASTA, writing a TSV to the given path.
# Usage: run_apex INPUT_FASTA OUTPUT_TSV
run_apex() {
    local input_fasta="$1"
    local output_tsv="$2"

    mkdir -p "$(dirname "${output_tsv}")"
    local prev_dir="${PWD}"
    cd "${APEX_DIR}"
    python predict.py "${input_fasta}" "${output_tsv}" all
    cd "${prev_dir}"
}

# Format a hydrophobicity float for use in directory names.
# e.g.  -0.75 → m0p75   0.0 → 0p0   0.25 → 0p25
fmt_hydro() {
    echo "$1" | tr '.' 'p' | tr '-' 'm'
}
