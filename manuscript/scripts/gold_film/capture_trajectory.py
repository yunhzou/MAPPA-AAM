"""Replay all branches of one saved gold cut, retaining atom-level events.

Input checkpoints must be trusted files produced by MAPPA, not downloads from
unknown sources. A 300-second watchdog bounds the replay process.
"""
from pathlib import Path
import argparse
import gzip
import json
import os
import signal
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))


def capture(archive, output, context=21):
    from mappa.artifacts import read_aam_checkpoint
    from mappa.search_trajectory import build_trajectory
    from mappa.viewers import growth_trace_html
    aam = read_aam_checkpoint(archive)
    terminals = [t for t in aam.graph.terminals if aam.graph.states[t].context == context]
    document = build_trajectory([dict(aam=aam, context=context, terminals=terminals)],
                                title='Gold rearrangement: actual atom-by-atom growth',
                                key_atoms=[4, 15, 17])
    output.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(document, separators=(',', ':')) + '\n').encode()
    (output / 'trajectory.json.gz').write_bytes(gzip.compress(encoded, mtime=0))
    (output / 'trajectory.html').write_text(growth_trace_html(document))
    print(f'Captured all {len(terminals)} terminal branches in context {context}.', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--context', type=int, default=21)
    p.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    a = p.parse_args()
    if a.worker:
        capture(a.archive, a.output, a.context)
    else:
        child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), str(a.archive.resolve()),
                                  str(a.output.resolve()), '--context', str(a.context), '--worker'],
                                 start_new_session=True,
                                 env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
        try:
            raise SystemExit(child.wait(timeout=300))
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGTERM)
            child.wait()
            raise SystemExit('Replay stopped by the 300-second watchdog.')
