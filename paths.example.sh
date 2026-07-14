#!/usr/bin/env bash
# Copy to paths.local.sh in the repo root and edit. paths.local.sh is gitignored.

export OMEGAMP_DIR="${HOME}/omegamp-inference"
export BATTLEAMP_DIR="${HOME}/battleamp-snakemake"

# Optional.
# export SIM_PYTHON="${HOME}/miniforge3/envs/omegamp/bin/python"
# export GPU=0

# Activation for the env that has snakemake + BattleAMP (used by scoring scripts).
export BATTLEAMP_ACTIVATE="${HOME}/.venvs/battleamp-snakemake/bin/activate"