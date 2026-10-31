#!/bin/bash
#SBATCH --account=nemotron_omni_vision
#SBATCH --partition=cpu_short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-4
#SBATCH --cpus-per-task=60
#SBATCH --mem=240G
#SBATCH --time=03:45:00
#SBATCH --signal=B:USR1@180
#SBATCH --job-name=mappa-enzymemap
set -euo pipefail
ROOT=/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa
SOURCE="$ROOT/MAPPA-AAM"
cd "$SOURCE"
export PYTHONPATH="$SOURCE/src:$SOURCE/bench"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
export MAPPA_BUNDLE_ROOTS="{\"/scratch/fsw/portfolios/nvr/projects/nvr_lpr_agentic/users/yunhengz/mappa\":\"$ROOT\",\"/lustre/fsw/portfolios/nvr/users/yunhengz/mappa\":\"$ROOT\"}"
BATCH="${MAPPA_BATCH_ROOT:-$ROOT/enzymemap_pool_20261030/dfw160}/$SLURM_ARRAY_TASK_ID"
exec "$ROOT/venv/bin/python" bench/enzymemap.py "$BATCH/index.json" "$BATCH/result" --workers 56 --memory-gib 4 --window 128 --run-seconds 12600 --resume "$@"
