#!/bin/bash
#SBATCH --job-name=mappa_decode_families
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-31%16
set -euo pipefail
run=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/current_decode_20261023_v1
review=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/failure_review_20261023
export PYTHONPATH="$run/source/src" MAPPA_NATIVE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$review/venv/bin/python" "$run/decode_shard.py" --result "$run/results/search_memory" --inputs "$review/MAPPA_failure_examples_20261023/02_AAM_memory_PMC11894599_10_gpt-6-luna" --output "$run/shards" --shard "$SLURM_ARRAY_TASK_ID" --shards 32
