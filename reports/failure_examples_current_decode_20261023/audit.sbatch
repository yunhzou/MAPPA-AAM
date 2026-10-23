#!/bin/bash
#SBATCH --job-name=mappa_decode_audit
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
run=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/current_decode_20261023_v1
review=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/failure_review_20261023
python="$review/venv/bin/python"
inputs="$review/MAPPA_failure_examples_20261023/02_AAM_memory_PMC11894599_10_gpt-6-luna"
export PYTHONPATH="$run/source/src" MAPPA_NATIVE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$python" "$run/merge_shards.py" --run "$run" --inputs "$inputs"
"$python" "$run/audit_decode.py" --result "$run/results/search_memory_partitioned" --inputs "$inputs"
"$python" "$run/check_membership.py" --result "$run/results/search_memory_partitioned" --inputs "$inputs"
