#!/bin/bash
#SBATCH --job-name=mappa-controlled-golden
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:55:00
#SBATCH --array=0-999%50
set -euo pipefail
cd /lustre/fsw/portfolios/nvr/users/yunhengz/Projects/coordinate_alignment/golden_controlled_20260915
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
AAM_PY="$PWD/../validation_20260912/.venv/bin/python"
"$AAM_PY" -c 'import run; assert run.read(run.ROOT/"pilot-proof.json")["status"]=="passed"'
aam_case_offset="${AAM_CASE_OFFSET:-0}"
aam_case="$((SLURM_ARRAY_TASK_ID + aam_case_offset))"
if [[ -f "results/$aam_case/done.json" ]]; then exit 0; fi
exec "$AAM_PY" run.py case "$aam_case"
