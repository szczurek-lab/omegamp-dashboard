#!/usr/bin/env bash
# Shared config for OmegAMP dashboard reproduce scripts.
# Run scripts from the dashboard repo root.

DASHBOARD_DIR="$(pwd)"

[ -f paths.local.sh ] && source paths.local.sh

OMEGAMP_DIR="${OMEGAMP_DIR:-${HOME}/omegamp-inference}"
BATTLEAMP_DIR="${BATTLEAMP_DIR:-${HOME}/battleamp-snakemake}"
BATTLEAMP_ACTIVATE="${BATTLEAMP_ACTIVATE:-${HOME}/.venvs/battleamp-snakemake/bin/activate}"

GENERATE_SCRIPT="${OMEGAMP_DIR}/project/scripts/inference/generate_samples.py"
PREDICT_SCRIPT="${OMEGAMP_DIR}/project/scripts/inference/predict_sequences.py"
CHECKPOINT="${OMEGAMP_DIR}/models/generative_model.ckpt"

POSITIVE_PROTOTYPES="${DASHBOARD_DIR}/data/prototypes/curated-AMPs_samples.fasta"
NEGATIVE_PROTOTYPES="${DASHBOARD_DIR}/data/prototypes/curated-Non-AMPs_samples.fasta"

# Overridable so a host whose data lives elsewhere (e.g. bury uses data_new/) can
# pin it in the gitignored paths.local.sh: `export DATA_DIR=.../data_new`.
DATA_DIR="${DATA_DIR:-${DASHBOARD_DIR}/data}"

OPTIMAL_POSITIVE_TAU=0.7
OPTIMAL_POSITIVE_SIGMA=1.0
OPTIMAL_NEGATIVE_TAU=0.3
OPTIMAL_NEGATIVE_SIGMA=1.0

SIM_PYTHON="${SIM_PYTHON:-python}"
GPU="${GPU:-0}"

# Repeat each FASTA entry N times, grouped by entry.
repeat_fasta() {
    awk -v n="$2" '
        function flush() { if (h) for (i=1;i<=n;i++) print h "\n" s }
        /^>/ { flush(); h=$0; s=""; next }
              { s = s $0 }
        END  { flush() }
    ' "$1" > "$3"
}