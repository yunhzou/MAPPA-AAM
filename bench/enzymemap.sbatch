#!/bin/bash
#SBATCH --account=nvr_lpr_agentic
#SBATCH --partition=cpu
#SBATCH --qos=cpu-short
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-3
#SBATCH --cpus-per-task=52
#SBATCH --mem=200G
#SBATCH --time=03:45:00
#SBATCH --signal=B:USR1@180
#SBATCH --job-name=mappa-enzymemap
set -euo pipefail
SOURCE=/scratch/fsw/portfolios/nvr/projects/nvr_lpr_agentic/users/yunhengz/mappa/MAPPA-AAM
PYTHON=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/rebenchmark_20261027/venv/bin/python
INDEX=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/rebenchmark_20261027/enzymemap_fresh_20261030/index.json
OUTPUT=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/enzymemap_pool_20261030
cd "$SOURCE"
export PYTHONPATH="$SOURCE/src:$SOURCE/bench"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
BATCH="$OUTPUT/bins100/$SLURM_ARRAY_TASK_ID"
exec "$PYTHON" bench/enzymemap.py "$BATCH/index.json" "$BATCH/result" --workers 48 --memory-gib 4 --window 128 --run-seconds 12600 --resume "$@"
