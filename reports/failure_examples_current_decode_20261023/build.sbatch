#!/bin/bash
#SBATCH --job-name=mappa_current_build
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -euo pipefail
run=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/current_decode_20261023_v1
python=/lustre/fsw/portfolios/nvr/users/yunhengz/mappa/failure_review_20261023/venv/bin/python
mkdir -p "$run/source"
tar xzf "$run/source.tar.gz" -C "$run/source"
cd "$run/source"
"$python" setup.py build_ext --inplace
"$python" native/build_engine.py
export PYTHONPATH="$run/source/src" MAPPA_NATIVE=1
"$python" - <<'PY'
from pathlib import Path
import hashlib,json
import mappa
from mappa.growth import native
root=Path.cwd();manifest=json.loads((root.parent/'source_sha256.json').read_text())
assert all(hashlib.sha256((root/p).read_bytes()).hexdigest()==v for p,v in manifest.items())
assert native.available()
assert not hasattr(mappa,'select_rp_mappings')
(root.parent/'build-validation.json').write_text(json.dumps({'source_files':len(manifest),'all_hashes_match':True,'mappa':mappa.__file__,'native':native.available(),'extensions':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'src/mappa').glob('*.so')}},indent=2))
PY
