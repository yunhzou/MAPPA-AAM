#!/usr/bin/env bash
#SBATCH --partition=cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=360G
#SBATCH --time=24:00:00

set -euo pipefail

: "${MAPPA_ROOT:?set MAPPA_ROOT to an immutable source snapshot}"
: "${MAPPA_PYTHON:?set MAPPA_PYTHON to the project interpreter}"
: "${TS_MANIFEST:?set TS_MANIFEST to one tier manifest}"
: "${TS_OUTPUT_ROOT:?set TS_OUTPUT_ROOT}"
: "${TS_WORK_ROOT:?set TS_WORK_ROOT}"
: "${TS_BENCHMARK_ROOT:?set TS_BENCHMARK_ROOT}"
: "${TS_REACTION_ROOT:?set TS_REACTION_ROOT to explicit resolved-reaction records}"

export PYTHONPATH="$MAPPA_ROOT/src"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

"$MAPPA_PYTHON" \
  "$MAPPA_ROOT/tools/run_typed_ts_case.py" \
  --manifest "$TS_MANIFEST" \
  --work-root "$TS_WORK_ROOT" \
  --benchmark-root "$TS_BENCHMARK_ROOT" \
  --output-root "$TS_OUTPUT_ROOT" \
  --reaction-root "$TS_REACTION_ROOT"
