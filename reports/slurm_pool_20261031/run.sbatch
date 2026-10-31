#!/bin/bash
#SBATCH --account=nvr_lpr_agentic
#SBATCH --partition=cpu
#SBATCH --qos=cpu-short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --job-name=mappa-pool-test
set -euo pipefail
SOURCE=/scratch/fsw/portfolios/nvr/projects/nvr_lpr_agentic/users/yunhengz/mappa/MAPPA-AAM
PYTHON=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/rebenchmark_20261027/venv/bin/python
OUTPUT=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/slurm_pool_case446_20261031
cd "$SOURCE"
export PYTHONPATH="$SOURCE/src:$SOURCE/tests"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
"$PYTHON" setup.py build_ext --inplace
"$PYTHON" src/mappa/mappa_core/native/build_engine.py
"$PYTHON" reports/slurm_pool_20261031/run_case.py "$OUTPUT"
