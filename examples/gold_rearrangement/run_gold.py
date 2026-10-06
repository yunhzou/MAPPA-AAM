"""Reproduce the endpoint-only 65-atom gold example and verify two oxygen fates.

Run from an installed MAPPA checkout. The default 300-second watchdog wraps
search, catalogue construction and the two selected-witness checks together.
No intermediate geometry, atom correspondence or pathway is a search input.
"""
from pathlib import Path
import argparse
import json
import os
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
DATA = Path(__file__).resolve().parent / 'data'


def run(output, workers=4):
    from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam
    from mappa.artifacts import write_aam_checkpoint
    from mappa.event_patterns import SignedEventIndex
    from mappa.family_query import query_path
    from mappa.growth.native import available

    raw = [json.loads((DATA / name).read_text()) for name in ('reactant.json', 'product.json')]
    problem = AAMProblem(*(MolecularEndpoint(d['elements'], d['coordinates'], d['wbo']) for d in raw),
                         name='Gold-catalyzed rearrangement: AuPPh3, 2 to 9')
    config = AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1.0, sweep_cuts=True)
    # search_aam is the published baseline; competition requires a separate opt-in API.
    backend = 'reused_native' if available() else 'reference'
    start = time.perf_counter()
    aam = search_aam(problem, config, workers=workers, execution=backend)
    search_seconds = time.perf_counter() - start
    catalogue = aam.final_catalogue()
    output.mkdir(parents=True, exist_ok=True)
    write_aam_checkpoint(aam, output / 'aam.checkpoint')

    # Reference oxygen assignments enter only after the independent search.
    # All indices below are zero-based. The two methyls at C6 are equivalent.
    oxygen_references = {'a': {4: 23, 15: 25, 17: 4},
                         'b_and_c': {4: 4, 15: 25, 17: 23}}
    r_core = set(raw[0]['organic_heavy_atoms'])
    p_core = set(raw[1]['organic_heavy_atoms'])
    carbons = sorted(i for i in r_core if problem.reactant.elements[i] == 'C')
    index = SignedEventIndex(problem, threshold=.5, metal_threshold=.3)
    selected = {}
    for name, reference in oxygen_references.items():
        candidates = []
        for fi, family in enumerate(catalogue.families):
            m = dict(family.mapping)
            if any(m[i] != j for i, j in reference.items()):
                continue
            if {m[i] for i in r_core} != p_core:
                continue
            if not all((problem.reactant.wbo[i, j] > .5) ==
                       (problem.product.wbo[m[i], m[j]] > .5)
                       for i in carbons for j in carbons if i < j):
                continue
            events = index.describe([m[i] for i in range(problem.source_atom_count)])
            candidates.append((events['total'], fi, events))
        if not candidates:
            raise RuntimeError(f'No carbon-framework-preserving representative recovered for {name}')
        _, fi, events = min(candidates)
        family = catalogue.families[fi]
        family.validate_representative(problem)
        status, witness = query_path(family.as_path(problem), problem, dict(family.mapping),
                                    source_atoms=range(problem.source_atom_count), timeout_ms=10000)
        if status != 'recovered':
            raise RuntimeError(f'Witness certification returned {status}')
        selected[name] = dict(family=family.to_record(), witness=witness, events=events,
                              oxygen_reference=reference,
                              carbon_framework_preserving_representatives=len(candidates))
    assert len({v['events']['id'] for v in selected.values()}) == 2
    summary = dict(atoms=problem.source_atom_count, config=vars(config), backend=backend,
                   workers=workers, search_seconds=search_seconds,
                   elapsed_through_selected_checks_seconds=time.perf_counter() - start,
                   final_fragment_branches=len(catalogue.branches),
                   saved_families=len(catalogue.families), saved_paths=catalogue.path_count,
                   incomplete_paths=catalogue.incomplete_paths,
                   scope='Two certified selected witnesses; not exhaustive event decoding.',
                   selected=selected)
    (output / 'recovery.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'selected'}, indent=2), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('gold-output'))
    parser.add_argument('--workers', type=int, default=min(4, max(1, (os.cpu_count() or 2) - 1)))
    parser.add_argument('--watchdog-seconds', type=int, default=300)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        run(args.output, args.workers)
    else:
        env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
               'VECLIB_MAXIMUM_THREADS': '1'}
        child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--worker',
                                  '--output', str(args.output.resolve()), '--workers', str(args.workers)],
                                 env=env, start_new_session=True)
        try:
            raise SystemExit(child.wait(timeout=args.watchdog_seconds))
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGTERM)
            child.wait()
            raise SystemExit('Stopped by the explicit watchdog; no result claimed.')
