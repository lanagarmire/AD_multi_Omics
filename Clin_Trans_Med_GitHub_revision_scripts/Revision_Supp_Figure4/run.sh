#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
python_bin="${PYTHON_BIN:-python3}"
"$python_bin" original_protocol_runner.py
"$python_bin" modality_baselines.py
"$python_bin" integrated_plot.py
