#!/bin/bash
#SBATCH --account=nemotron_omni_vision
#SBATCH --partition=cpu_short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --job-name=mappa-cache-pilot
set -euo pipefail
ROOT=/lustre/fsw/portfolios/nemotron/projects/nemotron_omni_vision/users/yunhengz/mappa
cd "$ROOT/MAPPA-AAM"
export PYTHONPATH="$PWD/src:$PWD/bench:$PWD/tests"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
export MAPPA_BUNDLE_ROOTS="{\"/scratch/fsw/portfolios/nvr/projects/nvr_lpr_agentic/users/yunhengz/mappa\":\"$ROOT\",\"/lustre/fsw/portfolios/nvr/users/yunhengz/mappa\":\"$ROOT\"}"
OUTPUT="$ROOT/cache-pilot-$SLURM_JOB_ID"
"$ROOT/venv/bin/python" src/mappa/mappa_core/native/build_engine.py > "$OUTPUT-build.log" 2>&1
"$ROOT/venv/bin/python" -m pytest -q tests/test_reaction_cache.py tests/test_fragment_repair.py > "$OUTPUT-tests.log" 2>&1
exec "$ROOT/venv/bin/python" bench/cache_pilot.py --output "$OUTPUT" --workers 8 \
  --index "$ROOT/enzymemap_pool_20261030/dfw160/0/index.json" \
  --index "$ROOT/enzymemap_pool_20261030/dfw160/1/index.json" \
  --index "$ROOT/enzymemap_pool_20261030/dfw160/2/index.json" \
  --index "$ROOT/enzymemap_pool_20261030/dfw160/3/index.json" \
  --index "$ROOT/enzymemap_pool_20261030/dfw160/4/index.json"
