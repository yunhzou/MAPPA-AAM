"""Count distinct heavy-atom choices at the saved detected event minimum.

This is external analysis of immutable decoded responses. Hydrogen assignments
are existential: a heavy mapping needs one valid completion, not one count per
hydrogen permutation. Chemical endpoint equivalence is the benchmark's existing
paired-graph certificate. No matching, decoding or minimum optimization is run.
"""
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import z3

from compare_slap import reference_certificate
from decoded_reference import audit_inputs
from golden import _publish
from mappa.mappa_core.core.scheduler.bundles import Blob
from mappa.mappa_core.core.symmetry.constraints import SymbolicActions
from mappa.mappa_core.utils.io.records import decoded_events_from_record


class ProjectedChoices:
    """Project the saved ordered operations onto mappings and their saved score."""
    def __init__(self, table, event_count, deadline):
        self.table, self.deadline = table, deadline
        self.solver = z3.Solver()
        self.encoder = SymbolicActions(self.solver, fresh=True)
        program = table.program
        self.prefix = self.encoder.act([self.encoder.constant(i) for i in range(program.size)],
                                      program.group.generators, 'saved_decoded_choices')
        self.values = self.mapping(program, self.prefix)
        constant = sum(int(np.asarray(row).sum()) for row in table.fixed_boundary.values())
        terms = []
        for group in table._term_groups():
            for term in group:
                self.check()
                response = term.response
                if response.constant is not None:
                    constant += int(bool(response.constant))
                    continue
                a, b = (self.values[i] for i in term.atoms)
                condition = self.pair_condition(a, b,
                    lambda x, y: bool(response.lookup(table.index.p[x, y])))
                terms.append(z3.If(condition, 1, 0))
        if table.sweep_cut is not None:
            constant -= int(table.boundary_events.source_cut_counted)
            a, b = (self.values[i] for i in table.sweep_cut)
            condition = self.pair_condition(a, b,
                lambda x, y: tuple(sorted((x, y))) in table.product_boundary_set)
            terms.append(-z3.If(condition, 1, 0))
        self.solver.add(constant + z3.Sum(terms) == event_count)
        self.event_count = event_count

    def check(self):
        if time.perf_counter() >= self.deadline:
            raise TimeoutError('choice-count allowance exhausted')

    def remap(self, value, transform):
        return self.encoder.lookup(value, {i: self.encoder.constant(transform(i)) for i in value[1]})

    def mapping(self, response, prefix):
        self.check()
        if hasattr(response, 'actions'):
            values = {a: self.encoder.constant(b) for a, b in response.base}
            for ((kind, data), prepared), (start, atoms, degree) in zip(response.actions, response.frames):
                local = [self.remap(value, lambda p: p-start) for value in prefix[start:start+degree]]
                if kind == 'relation':
                    selected = self.mapping(prepared, local)
                    transport = {b: selected[a] for a, b in data.witness}
                else:
                    transport = {a: self.remap(value, lambda p: atoms[p])
                                 for a, value in zip(atoms, local)}
                values = {a: self.encoder.lookup(value, transport) for a, value in values.items()}
            return values
        n = len(response.witness)
        targets = [self.remap(value, lambda p: response.fragment.target_atoms[response.witness[p]])
                   for value in (prefix if not response.two_sided else
                                 [self.remap(v, lambda p: p-n) for v in prefix[n:]])]
        if response.two_sided:
            # Source and target operations are joint permutations, never free pools.
            inverse = []
            for atom in range(n):
                possible = [i for i, value in enumerate(prefix[:n]) if atom in value[1]]
                if len(possible) == 1:
                    inverse.append(self.encoder.constant(possible[0]))
                else:
                    value = self.encoder.symbol(possible)
                    self.solver.add(z3.Or(*(z3.And(value[0] == i, prefix[i][0] == atom) for i in possible)))
                    inverse.append(value)
            targets = [self.encoder.lookup(value, dict(enumerate(targets))) for value in inverse]
        return dict(zip(response.fragment.source_atoms, targets))

    def pair_condition(self, left, right, predicate):
        rows = []
        for a in left[1]:
            self.check()
            targets = [b for b in right[1] if predicate(a, b)]
            if targets:
                rows.append(z3.And(left[0] == a, z3.Or(*(right[0] == b for b in targets))))
        return z3.simplify(z3.Or(*rows))

    def next(self, heavy):
        self.check()
        self.solver.set(timeout=max(1, int(1000*(self.deadline-time.perf_counter()))))
        status = self.solver.check()
        if status == z3.unknown:
            raise TimeoutError(self.solver.reason_unknown())
        if status == z3.unsat:
            return None
        model = self.solver.model()
        def number(value):
            return value if isinstance(value, int) else model.eval(value, model_completion=True).as_long()
        prefix = np.asarray([number(value[0]) for value in self.prefix])
        mapping = self.table.program.mapping(prefix)
        if mapping != {a: number(value[0]) for a, value in self.values.items()}:
            raise AssertionError('Symbolic mapping differs from saved ordered action replay')
        vector = self.table.index.vector(mapping)
        if (self.table.program.fixed_events + self.table.program.score(prefix, self.table.program.pairs)
                != self.event_count or self.table.evaluate(vector)['total'] != self.event_count):
            raise AssertionError('Choice does not have the saved requested event count')
        self.table.index._vectors([vector])
        pairs = tuple(sorted((a, b) for a, b in mapping.items() if a in heavy))
        self.solver.add(z3.Or(*(self.values[a][0] != b for a, b in pairs)))
        return pairs


def count_choices(decoded_directions, reference, *, seconds=60., max_mappings=100000):
    """Exact union at the detected score, or an explicit lower bound on choices."""
    started = time.perf_counter()
    deadline = started + seconds
    minimum = min(decoded.minimum_events for _, decoded in decoded_directions)
    expected = reference_certificate(reference['features'], reference['mapping'])
    classes, labelled = {}, set()
    contexts, completed, queries = 0, 0, 0
    exhausted = None
    try:
        for reverse, decoded in decoded_directions:
            heavy = set(reference['features'][int(reverse)]['heavy'])
            seen_tables = set()
            for context in decoded.window.contexts:
                if decoded.window.bounds[context] > minimum:
                    continue
                key = decoded._table_key(*context)
                if key in seen_tables:
                    continue
                seen_tables.add(key)
                contexts += 1
                if time.perf_counter() >= deadline:
                    raise TimeoutError('choice-count allowance exhausted')
                table = decoded._prepare_table(*context[:2], key)
                projected = ProjectedChoices(table, minimum, deadline)
                while True:
                    pairs = projected.next(heavy)
                    queries += 1
                    if pairs is None:
                        break
                    pairs = tuple(sorted((b, a) for a, b in pairs)) if reverse else pairs
                    if pairs not in labelled:
                        labelled.add(pairs)
                        certificate = reference_certificate(reference['features'], pairs)
                        classes.setdefault(certificate, pairs)
                    if queries >= max_mappings:
                        raise TimeoutError('mapping-count allowance exhausted')
                completed += 1
                del projected, table
    except TimeoutError as error:
        exhausted = str(error)
    complete = exhausted is None
    recovered = expected in classes
    if complete and not classes:
        raise AssertionError('Saved detected minimum has no replayable mapping')
    return dict(status='complete' if complete else 'incomplete', detected_minimum=minimum,
        saved_minimum_proven=all(d.minimum_proven or d.minimum_lower_bound >= minimum
                                for _, d in decoded_directions),
        choices=len(classes) if complete else None, choices_lower_bound=len(classes),
        reference_at_minimum=recovered if complete or recovered else None,
        random_success=(int(recovered)/len(classes) if complete else None),
        labelled_heavy_mappings=len(labelled), contexts_visited=contexts,
        contexts_completed=completed, solver_queries=queries, reason=exhausted,
        seconds=time.perf_counter()-started, representatives=list(classes.values()))


def analyze(payload, workspace):
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != payload['runner_sha256']:
        raise ValueError('Analysis source changed after planning')
    decoded = []
    for row in payload['directions']:
        saved = workspace.load(Blob(**row['decoded']))
        decoded.append((row['reversed'], decoded_events_from_record(saved.record)))
    reference = workspace.load(Blob(**payload['directions'][0]['reference']))
    return count_choices(decoded, reference, seconds=payload['seconds'], max_mappings=payload['max_mappings'])


def run(args):
    from mappa.execution import BundleStore, WorkPool
    from mappa.execution.identity import code_identity
    from mappa.mappa_core.core.scheduler import Task
    rows, sources, _ = audit_inputs(args.campaign)
    grouped = defaultdict(list)
    for row in rows:
        if args.cases is None or row['case'] in args.cases:
            grouped[row['case']].append(row)
    if not grouped or any(sorted(r['reversed'] for r in pair) != [False, True] for pair in grouped.values()):
        raise ValueError('Require exactly both directions of every requested reaction')
    args.output.mkdir(parents=True, exist_ok=False)
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report = dict(schema='mappa.minimum_choices/v1', original_hashes=sources,
        code_identity=code_identity(), runner_sha256=digest, expected=len(grouped),
        policy='Distinct paired heavy-atom chemical mapping classes at the lowest saved detected model count across both directions',
        hydrogen_policy='One full valid completion suffices; hydrogen permutations do not multiply choices',
        resources=dict(seconds=args.seconds, task_seconds=args.task_seconds,
                       max_mappings=args.max_mappings, workers=args.workers, memory_gib=args.memory_gib), results=[])
    store = BundleStore(args.output/'bundles')
    factory = WorkPool.slurm if args.slurm else WorkPool.local
    started = time.perf_counter()
    with factory(workers=args.workers, store=store, memory_bytes=int(args.memory_gib*2**30),
                 case_seconds=args.task_seconds, shared_cache=False) as pool:
        report['ready_workers'] = len(pool.queue.ready_workers())
        waiting = iter(sorted(grouped.items()))
        pending, done = {}, False
        publish_at = 0.
        while pending or not done:
            from itertools import islice
            admitted = list(islice(waiting, 128-len(pending)))
            if len(pending) < 128:
                done = not admitted
            for case, directions in admitted:
                key = str(case)
                pending[key] = case
                pool.queue.add([Task(key, 'minimum_choices:analyze', dict(directions=directions,
                    seconds=args.seconds, max_mappings=args.max_mappings, runner_sha256=digest),
                    case=key, case_seconds=args.task_seconds, memory_bytes=int(args.memory_gib*2**30))])
            pool.check_workers()
            for key, state in pool.queue.results(tuple(pending)).items():
                if state['state'] in ('pending', 'running'):
                    continue
                row = dict(case=pending.pop(key), state=state['state'])
                measured = state['result'] or {}
                row['metrics'] = {k:v for k,v in measured.items() if k.endswith(('_seconds','_bytes'))}
                if state['state'] == 'complete':
                    value = store.load(measured['blob'])
                    row.update({k:v for k,v in value.items() if k != 'representatives'})
                    row['evidence'] = asdict(measured['blob'])
                else:
                    row.update(status='incomplete', error=measured.get('error'))
                report['results'].append(row)
                if state['state'] == 'error':
                    pool.queue.finish_input(cancel_pending=True)
                    _publish(args.output/'results.json', report)
                    raise RuntimeError('Choice-count analysis failed: '+str(row))
            if time.perf_counter() >= publish_at or not pending:
                report.update(elapsed_seconds=time.perf_counter()-started,
                    status_counts=dict(Counter(r['status'] for r in report['results'])))
                _publish(args.output/'results.json', report)
                print(json.dumps(dict(completed=len(report['results']), expected=report['expected'],
                                      status_counts=report['status_counts'])), flush=True)
                pool.checkpoint()
                publish_at = time.perf_counter()+30
            if pending:
                time.sleep(.1)
    report['source_unchanged'] = all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p,h in sources.items())
    if not report['source_unchanged'] or len(report['results']) != report['expected']:
        raise AssertionError('Original source or complete analysis scope changed')
    report['elapsed_seconds'] = time.perf_counter()-started
    _publish(args.output/'results.json', report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cases', type=int, nargs='+')
    parser.add_argument('--seconds', type=float, default=60.)
    parser.add_argument('--task-seconds', type=float, default=180.)
    parser.add_argument('--max-mappings', type=int, default=100000)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--memory-gib', type=float, default=4.)
    parser.add_argument('--slurm', action='store_true')
    run(parser.parse_args())
