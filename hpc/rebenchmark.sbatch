#!/bin/bash
#SBATCH --job-name=mappa_golden
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=03:50:00
set -euo pipefail
repo=$1
run=$2
export PYTHONPATH="$repo/src:$repo/bench"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
"$run/../venv/bin/python" -m rebenchmark.campaign drain --run "$run"
