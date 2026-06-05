#!/usr/bin/env bash
# mapseq2-py installer for Linux / WSL2 (Ubuntu).
#
# Usage:
#   bash install.sh             # install into a fresh `mapseq` conda env
#   ENV_NAME=foo bash install.sh # install into a custom-named conda env
#
# What it does:
#   1. Verifies a working conda is on PATH (installs miniforge if not).
#   2. Creates conda env `mapseq` from environment.yml
#      (bowtie2 + pigz + scientific Python + this package in editable mode).
#   3. Smoke-tests that `mapseq2 --help` and `bowtie2 --version` work.
#
# Requirements assumed:
#   - bash, curl, basic build tools (usually present on Ubuntu/WSL by default)
#   - internet access for conda + bioconda package downloads

set -euo pipefail

ENV_NAME="${ENV_NAME:-mapseq}"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

log() { printf '\033[1;34m[install.sh]\033[0m %s\n' "$*"; }
err() { printf '\033[1;31m[install.sh ERROR]\033[0m %s\n' "$*" >&2; }

# ---------- 1. Conda ---------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
  log "conda not found on PATH; installing Miniforge to \$HOME/miniforge3"
  TMP_INSTALLER="$(mktemp --suffix=.sh)"
  curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -o "$TMP_INSTALLER"
  bash "$TMP_INSTALLER" -b -p "$HOME/miniforge3"
  rm -f "$TMP_INSTALLER"
  # shellcheck disable=SC1091
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
  conda init bash >/dev/null
  log "Miniforge installed. Restart your shell after this script finishes."
else
  CONDA_BASE="$(conda info --base)"
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

# ---------- 2. Environment ---------------------------------------------------
if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  log "conda env '$ENV_NAME' already exists; updating from environment.yml"
  conda env update -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml" --prune
else
  log "creating conda env '$ENV_NAME' from environment.yml"
  conda env create -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml"
fi

conda activate "$ENV_NAME"

# ---------- 3. Smoke test ----------------------------------------------------
log "verifying installation"
bowtie2 --version | head -1
pigz --version 2>&1 | head -1
python -c "import mapseq2; print('mapseq2 import: OK')"
mapseq2 --help >/dev/null && log "mapseq2 CLI: OK"

log "DONE."
cat <<EOF

Next steps:
  1. Activate the env in new shells:
       conda activate $ENV_NAME

  2. Try the example config:
       cp configs/example_config.yaml my_config.yaml
       cp configs/example_samplesheet.csv my_samples.csv
       # edit my_config.yaml + my_samples.csv to point at your FASTQs
       mapseq2 -v run --config my_config.yaml

  3. See MANUAL.md for the full workflow.
EOF
