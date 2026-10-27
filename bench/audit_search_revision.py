"""Audit one directional search from an explicitly selected source checkout.

Historical checkouts are diagnostic comparisons, not supported runtime modes.
Optional phase timers measure exclusive process CPU; no matching is changed.
"""
import argparse
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import asdict
import hashlib
import gc
import json
from pathlib import Path
import resource
import sys
import time


class Phases:
    def __init__(self):
        self.stack = []
        self.cpu = defaultdict(float)
        self.calls = defaultdict(int)

    @contextmanager
    def measure(self, name):
        frame = [time.process_time(), 0.]
        self.stack.append(frame)
        try:
            yield
        finally:
            elapsed = time.process_time() - frame[0]
            self.stack.pop()
            self.cpu[name] += elapsed - frame[1]
            self.calls[name] += 1
            if self.stack:
                self.stack[-1][1] += elapsed

    def wrap(self, owner, name, category, *, class_method=False):
        original = getattr(owner, name)
        if class_method:
            def wrapped(_cls, *args, **kwargs):
                with self.measure(category):
                    return original(*args, **kwargs)
            setattr(owner, name, classmethod(wrapped))
        else:
            def wrapped(*args, **kwargs):
                with self.measure(category):
                    return original(*args, **kwargs)
            setattr(owner, name, wrapped)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reverse', action='store_true')
    parser.add_argument('--phases', action='store_true')
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo / 'src'))
    from mappa import AAMProblem, AAMSearchConfig, MolecularEndpoint
    from mappa import aam, artifacts, _engine
    from mappa.search_graph import AAMSearchGraph
    from mappa.conditioned_symmetry import ConditionedSymmetryWorkspace
    from mappa.matcher.canonical import _CandidateAutomorphismCanonicalizer

    args.output.mkdir(parents=True, exist_ok=False)
    raw = json.loads(args.input.read_text())
    endpoints = [MolecularEndpoint(**raw[k]) for k in ('reactant', 'product')]
    if args.reverse:
        endpoints.reverse()
    problem = AAMProblem(*endpoints, raw['name'])
    config = AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1., random_seed=42)
    phases = Phases()
    rows, workspaces = [], []
    original_append = aam._GrowthCounts.append
    def append(self, row):
        rows.append(row)
        original_append(self, row)
    aam._GrowthCounts.append = append
    original_init = ConditionedSymmetryWorkspace.__init__
    def initialize(self, *a, **kw):
        original_init(self, *a, **kw)
        workspaces.append(self)
    ConditionedSymmetryWorkspace.__init__ = initialize
    if args.phases:
        for owner, name, category, class_method in (
            (aam, '_initialize_search', 'initialization', False),
            (aam, '_search_cut', 'search_other', False),
            (_engine, 'search_fragments', 'native_search_and_export', False),
            (AAMSearchGraph, 'from_record', 'python_graph_construction', True),
            (AAMSearchGraph, 'combine', 'graph_combination', True),
            (artifacts, 'write_raw_cut', 'raw_checkpoint_writes', False),
            (artifacts, 'write_aam_checkpoint', 'final_checkpoint_write', False),
            (aam, 'finalize_graph_symmetry', 'symmetry_bookkeeping', False),
            (ConditionedSymmetryWorkspace, 'conditioned_generators', 'conditioned_groups', False),
        ):
            phases.wrap(owner, name, category, class_method=class_method)
        if hasattr(_CandidateAutomorphismCanonicalizer, '_growth_roles'):
            phases.wrap(_CandidateAutomorphismCanonicalizer, '_growth_roles', 'continuation_roles')
        collection = []
        def track_gc(event, info):
            if event == 'start':
                context = phases.measure('cyclic_gc')
                context.__enter__()
                collection.append(context)
            else:
                collection.pop().__exit__(None, None, None)
        gc.callbacks.append(track_gc)

    cpu, wall = time.process_time(), time.perf_counter()
    with phases.measure('other'):
        result = aam.search_aam(problem, config, workers=1, execution='reused_native',
                                intermediate_dir=args.output / 'cuts', archive_format='checkpoint')
    cpu, wall = time.process_time() - cpu, time.perf_counter() - wall
    if args.phases:
        gc.callbacks.remove(track_gc)
    graph = result.graph
    source_hashes = {str(p.relative_to(args.repo)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for folder in ('src/mappa', 'native/src')
                     for p in (args.repo / folder).rglob('*')
                     if p.is_file() and p.suffix in ('.py', '.h', '.cpp', '.so')}
    native_counts = {key: sum(r.get(key, 0) for r in rows)
                     for key in ('growth_calls', 'native_schedule_and_export_cpu', 'python_graph_cpu')}
    record = dict(config=asdict(config), direction='P_to_R' if args.reverse else 'R_to_P',
                  cpu_seconds=cpu, wall_seconds=wall, phase_instrumentation=args.phases,
                  phase_cpu=dict(phases.cpu), phase_calls=dict(phases.calls),
                  native_counts=native_counts, repair_stats=aam._SEARCH_REPAIR.stats(),
                  target_certificate_stats=list(aam._SEARCH_REPAIR.target_view.graph.memo_stats()),
                  symmetry_stats=[w.stats() for w in workspaces], metrics=asdict(result.metrics),
                  states=len(graph.states), transitions=len(graph.transitions),
                  terminals=len(graph.terminals), capped_stops=sum(s.reason == 'capped' for s in graph.stops),
                  graph_sha256=hashlib.sha256(graph.to_json().encode()).hexdigest(),
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  source_sha256=source_hashes,
                  peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                                 * (1 if sys.platform == 'darwin' else 1024))
    (args.output / 'measurement.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k:record[k] for k in ('cpu_seconds', 'phase_cpu', 'native_counts',
                                          'states', 'transitions', 'terminals')}), flush=True)


if __name__ == '__main__':
    main()
