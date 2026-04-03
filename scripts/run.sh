#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${T2LEAD_CONDA_ENV:-t2lead}"
CONDA_ROOT="${CONDA_PREFIX:-/root/miniconda3}"
ENV_BIN="${CONDA_ROOT}/envs/${CONDA_ENV}/bin"

if [ ! -d "$ENV_BIN" ]; then
    echo "ERROR: conda env '${CONDA_ENV}' not found at ${ENV_BIN}" >&2
    exit 1
fi

export PATH="${ENV_BIN}:${PATH}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export PYTHONPATH="${PROJECT_ROOT}/src${PYTHONPATH:+:$PYTHONPATH}"

exec "${ENV_BIN}/python" -m drugpipe.pipeline "$@"
