#!/bin/bash
#SBATCH --job-name=mappa_current_decode
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-1
set -euo pipefail
run=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/current_decode_20261023_v1
review=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/failure_review_20261023
export PYTHONPATH="$run/source/src" MAPPA_NATIVE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cases=(rp_timeout search_memory)
case_name=${cases[$SLURM_ARRAY_TASK_ID]}
"$review/venv/bin/python" "$run/run_decode.py" --case "$case_name" --inputs "$review/MAPPA_failure_examples_20261023" --output "$run/results/$case_name" --memory-gib 5
