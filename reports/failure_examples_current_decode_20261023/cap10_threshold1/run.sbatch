#!/bin/bash
#SBATCH --job-name=mappa_cap10_thr1
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
run=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/current_decode_20261023_v1
review=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/failure_review_20261023
export PYTHONPATH="$run/source/src" MAPPA_NATIVE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
experiment="$run/cap10_threshold1"
"$review/venv/bin/python" "$experiment/run_decode.py" --case search_memory --inputs "$review/MAPPA_failure_examples_20261023" --output "$experiment/result" --memory-gib 5
