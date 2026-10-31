#!/bin/bash
#SBATCH --account=nvr_lpr_agentic
#SBATCH --partition=cpu
#SBATCH --qos=cpu-short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=96G
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
# Verify the prepared-data adapter through every stage before opening the full pool.
"$PYTHON" bench/enzymemap.py "$INDEX" "$OUTPUT/input-check" --local --workers 2 --window 4 --limit 2
exec "$PYTHON" bench/enzymemap.py "$INDEX" "$OUTPUT/full" --workers 32 --window 128 --run-seconds 12600
