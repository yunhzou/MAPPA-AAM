#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p build
manuscript_python="${MANUSCRIPT_PYTHON:-python3}"
"$manuscript_python" scripts/build_figures.py
tectonic --keep-logs --keep-intermediates --outdir build main.tex
cp build/main.pdf manuscript.pdf
"$manuscript_python" scripts/validate_outputs.py
"$manuscript_python" scripts/package_bundle.py
