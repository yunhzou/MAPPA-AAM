#!/bin/bash
#SBATCH --account=nemotron_omni_vision
#SBATCH --partition=cpu_short
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=54G
#SBATCH --array=0-4
#SBATCH --time=01:00:00
#SBATCH --signal=B:USR1@180
#SBATCH --job-name=mappa-golden
set -euo pipefail
MAPPA_ROOT=/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa
cd "$MAPPA_ROOT/MAPPA-AAM"
export PYTHONPATH="$PWD/src:$PWD/bench"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
: "${MAPPA_GOLDEN_ROOT:?Set MAPPA_GOLDEN_ROOT to the prepared campaign directory}"
export MAPPA_BUNDLE_ROOTS="{\"${MAPPA_GOLDEN_INPUT_ROOT:-$MAPPA_GOLDEN_ROOT}\":\"$MAPPA_GOLDEN_ROOT\"}"
BATCH="$MAPPA_GOLDEN_ROOT/$SLURM_ARRAY_TASK_ID"
exec "$MAPPA_ROOT/venv/bin/python" bench/enzymemap.py "$BATCH/index.json" "$BATCH/result" \
  --workers 25 --memory-gib 2 --window 128 --run-seconds 3000 "$@"
