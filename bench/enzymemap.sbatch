#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --signal=B:USR1@180
#SBATCH --job-name=mappa-benchmark
set -euo pipefail

: "${MAPPA_SOURCE:?Set the checkout directory}"
: "${MAPPA_PYTHON:?Set the Python interpreter path}"
: "${MAPPA_BATCH_ROOT:?Set the prepared campaign directory}"
: "${MAPPA_EXPECTED_CODE_IDENTITY:?Set the validated source and native engine identity}"
: "${MAPPA_WORKERS_PER_NODE:?Set the number of compute workers per node}"
: "${SLURM_ARRAY_TASK_ID:?Submit this launcher as a Slurm job array}"

cd "$MAPPA_SOURCE"
export PYTHONPATH="$MAPPA_SOURCE/src:$MAPPA_SOURCE/bench"
export MAPPA_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export TMPDIR="${SLURM_TMPDIR:-/tmp}"

"$MAPPA_PYTHON" - <<'IDENTITY'
import os
from mappa.execution.identity import code_identity

if code_identity() != os.environ['MAPPA_EXPECTED_CODE_IDENTITY']:
    raise SystemExit('Worker source or native engine changed')
IDENTITY

BATCH="$MAPPA_BATCH_ROOT/$SLURM_ARRAY_TASK_ID"
exec "$MAPPA_PYTHON" bench/enzymemap.py "$BATCH/index.json" "$BATCH/result" \
  --workers "$MAPPA_WORKERS_PER_NODE" "$@"
