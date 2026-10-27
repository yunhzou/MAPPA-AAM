"""Measure one specified directional search, including saved checkpoints.

Run individual cases sequentially in fresh processes. This is not a dataset
campaign launcher. Imports and result hashing are outside the CPU measurement.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

from mappa import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_aam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--reverse', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    raw = json.loads(args.input.read_text())
    endpoints = [MolecularEndpoint(**raw[k]) for k in ('reactant', 'product')]
    if args.reverse:
        endpoints.reverse()
    problem = AAMProblem(*endpoints, raw['name'])
    config = AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1.,
                             random_seed=42)
    root = Path(__file__).resolve().parents[1]
    sources = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
               for folder in ('src/mappa', 'native/src')
               for p in (root / folder).rglob('*')
               if p.is_file() and p.suffix in ('.py', '.cpp', '.h', '.so')}
    cpu, wall = time.process_time(), time.perf_counter()
    result = search_aam(problem, config, workers=1, execution='reused_native',
                        intermediate_dir=args.output / 'cuts', archive_format='checkpoint')
    cpu, wall = time.process_time() - cpu, time.perf_counter() - wall
    graph = result.graph.to_json().encode()
    record = dict(input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  direction='P_to_R' if args.reverse else 'R_to_P', config=asdict(config),
                  cpu_seconds=cpu, wall_seconds=wall,
                  peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                                 * (1 if sys.platform == 'darwin' else 1024),
                  graph_sha256=hashlib.sha256(graph).hexdigest(),
                  terminals=len(result.graph.terminals), source_sha256=sources,
                  base_commit=subprocess.check_output(['git', '-C', str(root),
                                                       'rev-parse', 'HEAD'], text=True).strip())
    (args.output / 'measurement.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: v for k, v in record.items() if k != 'source_sha256'}), flush=True)


if __name__ == '__main__':
    main()
