#!/bin/bash
#SBATCH --job-name=mappa-small-pilot
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
cd /lustre/fsw/portfolios/nvr/users/yunhengz/Projects/coordinate_alignment/golden_controlled_20260915
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
AAM_PY="$PWD/../validation_20260912/.venv/bin/python"
"$AAM_PY" small_pilot.py
