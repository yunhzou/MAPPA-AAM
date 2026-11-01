#!/bin/bash
#SBATCH --account=nemotron_omni_vision
#SBATCH --partition=cpu_short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-3
#SBATCH --cpus-per-task=52
#SBATCH --mem=104G
#SBATCH --time=03:45:00
#SBATCH --signal=B:USR1@180
#SBATCH --job-name=mappa-enzymemap
set -euo pipefail
ROOT=/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa
SOURCE="$ROOT/MAPPA-AAM"
cd "$SOURCE"
export PYTHONPATH="$SOURCE/src:$SOURCE/bench"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export MAPPA_NATIVE=1
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
export MAPPA_BUNDLE_ROOTS="{\"/scratch/fsw/portfolios/nvr/projects/nvr_lpr_agentic/users/yunhengz/mappa\":\"$ROOT\",\"/lustre/fsw/portfolios/nvr/users/yunhengz/mappa\":\"$ROOT\"}"
: "${MAPPA_BATCH_ROOT:?Set MAPPA_BATCH_ROOT to the prepared campaign directory}"
: "${MAPPA_EXPECTED_CODE_IDENTITY:?Set the validated source and native engine identity}"
"$ROOT/venv/bin/python" -c 'import os; from mappa.execution.identity import code_identity; assert code_identity() == os.environ["MAPPA_EXPECTED_CODE_IDENTITY"], "Worker source or native engine changed"'
BATCH="$MAPPA_BATCH_ROOT/$SLURM_ARRAY_TASK_ID"
exec "$ROOT/venv/bin/python" bench/enzymemap.py "$BATCH/index.json" "$BATCH/result" \
  --workers "${MAPPA_WORKERS_PER_NODE:-50}" --memory-gib 2 \
  --window "${MAPPA_ACTIVE_DIRECTIONS:-128}" --run-seconds 9600 "$@"
