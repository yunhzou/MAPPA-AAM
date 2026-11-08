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
    def __init__(self, table, event_count, deadline, *, heavy=None, known=()):
        self.table, self.deadline = table, deadline
        self.solver = z3.Solver()
        self.encoder = SymbolicActions(self.solver, fresh=True)
        program = table.program
        self.prefix = self.encoder.act([self.encoder.constant(i) for i in range(program.size)],
                                      program.group.generators, 'saved_decoded_choices')
        self.values = self.mapping(program, self.prefix)
        self.domain_solver = z3.Solver()
        self.domain_values = {a: z3.FreshInt('heavy_image') for a in self.values
                              if heavy is None or a in heavy}
        for a,value in self.domain_values.items():
            self.domain_solver.add(z3.Or(*(value == p for p in self.values[a][1])))
        self.domain_solver.add(z3.Distinct(*self.domain_values.values()))
        for pairs in known:
            if set(dict(pairs)) == set(self.domain_values) and all(b in self.values[a][1] for a,b in pairs):
                self.exclude(pairs)
        self.event_count = event_count
        self.covered = self.domain_covered()
        if self.covered:
            return
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

    def exclude(self, pairs):
        self.solver.add(z3.Or(*(self.values[a][0] != b for a,b in pairs)))
        self.domain_solver.add(z3.Or(*(self.domain_values[a] != b for a,b in pairs)))

    def domain_covered(self):
        """A relaxed injective image domain can certify no unseen heavy output.

        This ignores event and symmetry restrictions, so only UNSAT can skip the
        actual query. It avoids proving absence through irrelevant H operations.
        """
        self.check()
        self.domain_solver.set(timeout=max(1,min(50,int(1000*(self.deadline-time.perf_counter())))))
        return self.domain_solver.check() == z3.unsat

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
        if self.covered or self.domain_covered():
            return None
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
        self.exclude(pairs)
        return pairs


def chemical_orbit(pairs, generators, *, limit=4096, deadline=float('inf')):
    """Bounded heavy-only images of a class already witnessed at the minimum.

    These images are exclusions, not extra predictions: every image has the same
    paired chemical graph. Stopping early only leaves more duplicate work.
    """
    seen, pending = {pairs}, [pairs]
    while pending and len(seen) < limit and time.perf_counter() < deadline:
        old = pending.pop()
        for side, operations in enumerate(generators):
            for g in operations:
                image = tuple(sorted((a,g[b]) if side else (g[a],b) for a,b in old))
                if image not in seen:
                    seen.add(image); pending.append(image)
                    if len(seen) >= limit:
                        return seen
    return seen


def heavy_images(program, heavy, *, deadline, limit=4096):
    """A small exact heavy projection of ordered actions, independent of events.

    Use only as a coverage certificate for classes already found at the minimum.
    Large projections are left to the existing bounded event query.
    """
    initial = dict(program.base)
    atoms = tuple(a for a in initial if a in heavy)
    states = {tuple(initial[a] for a in atoms)}
    for (kind,data),prepared in program.actions:
        if time.perf_counter() >= deadline:
            return None
        if kind == 'relation':
            selected = prepared.mapping(np.arange(prepared.size))
            transport = {b:selected[a] for a,b in data.witness}
            states = {tuple(transport.get(p,p) for p in state) for state in states}
            fragment = prepared.fragment
            groups = []
            if prepared.two_sided:
                source = fragment.source_atoms
                groups.append([{selected[a]:selected[source[g[i]]] for i,a in enumerate(source)
                                if i != g[i]} for g in fragment.source_symmetry.generators])
            target = fragment.target_atoms
            groups.append([{a:target[g[i]] for i,a in enumerate(target) if i != g[i]}
                           for g in fragment.target_symmetry.generators])
        else:
            groups = [[{a:b for a,b in enumerate(g) if a != b} for g in data]]
        for generators in groups:
            pending = list(states)
            while pending:
                if time.perf_counter() >= deadline:
                    return None
                old = pending.pop()
                for g in generators:
                    new = tuple(g.get(p,p) for p in old)
                    if new not in states:
                        if len(states) >= limit:
                            return None
                        states.add(new); pending.append(new)
    return {tuple(zip(atoms,state)) for state in states}


def count_choices(decoded_directions, reference, *, seconds=60., max_mappings=100000):
    """Exact union at the detected score, or an explicit lower bound on choices."""
    started = time.perf_counter()
    deadline = started + seconds
    minimum = min(decoded.minimum_events for _, decoded in decoded_directions)
    expected = reference_certificate(reference['features'], reference['mapping'])
    from golden import endpoint_generators
    generators = reference.get('generators')
    if generators is None:
        generators = [endpoint_generators(f,1+max(f['heavy'],default=-1)) for f in reference['features']]
    classes, labelled, excluded = {}, set(), set()
    contexts, completed, queries = 0, 0, 0
    exhausted = None
    try:
        for reverse, decoded in decoded_directions:
            heavy = set(reference['features'][int(reverse)]['heavy'])
            seen_tables, covered_programs = set(), set()
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
                if key[:2] in covered_programs:
                    completed += 1
                    continue
                table = decoded._prepare_table(*context[:2], key)
                if excluded:
                    images = heavy_images(table.program,heavy,deadline=min(deadline,time.perf_counter()+.05))
                    if images is not None:
                        images = {tuple(sorted((b,a) for a,b in pairs)) if reverse else tuple(sorted(pairs))
                                  for pairs in images}
                        if images <= excluded:
                            covered_programs.add(key[:2])
                            completed += 1
                            continue
                known = [tuple(sorted((b,a) for a,b in pairs)) if reverse else pairs for pairs in excluded]
                projected = ProjectedChoices(table, minimum, deadline, heavy=heavy, known=known)
                while True:
                    pairs = projected.next(heavy)
                    queries += 1
                    if pairs is None:
                        break
                    pairs = tuple(sorted((b, a) for a, b in pairs)) if reverse else pairs
                    if pairs not in labelled:
                        labelled.add(pairs)
                        certificate = reference_certificate(reference['features'], pairs)
                        if certificate not in classes:
                            classes[certificate] = pairs
                            images = chemical_orbit(pairs,generators,deadline=deadline)
                            excluded.update(images)
                            for image in images:
                                oriented = tuple(sorted((b,a) for a,b in image)) if reverse else image
                                if (set(dict(oriented)) == set(projected.domain_values)
                                        and all(b in projected.values[a][1] for a,b in oriented)):
                                    projected.exclude(oriented)
                        excluded.add(pairs)
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
        labelled_heavy_mappings=len(labelled), equivalent_labels_excluded=len(excluded), contexts_visited=contexts,
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
    rows, sources, unmatched_heavy = audit_inputs(args.campaign)
    if not unmatched_heavy:
        raise ValueError('This metric requires exact heavy-atom references, including unmatched atoms')
    grouped = defaultdict(list)
    for row in rows:
        if args.cases is None or row['case'] in args.cases:
            grouped[row['case']].append(row)
    if not grouped or any(sorted(r['reversed'] for r in pair) != [False, True] for pair in grouped.values()):
        raise ValueError('Require exactly both directions of every requested reaction')
    previous, completed = None, []
    if args.remaining_from is not None:
        previous = json.loads(args.remaining_from.read_bytes())
        if any(previous['original_hashes'].get(p) != h for p,h in sources.items()):
            raise ValueError('Continuation inputs differ from the original saved analysis')
        if (previous['schema'] != 'mappa.minimum_choices/v1' or previous['expected'] != len(grouped)
                or not {r['case'] for r in previous['results']} <= set(grouped)
                or len({r['case'] for r in previous['results']}) != len(previous['results'])):
            raise ValueError('Continuation must belong to the same reaction scope')
        completed = [r for r in previous['results'] if r['status'] == 'complete']
        sources[str(args.remaining_from.resolve())] = hashlib.sha256(args.remaining_from.read_bytes()).hexdigest()
    args.output.mkdir(parents=True, exist_ok=False)
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report = dict(schema='mappa.minimum_choices/v1', original_hashes=sources,
        code_identity=code_identity(), runner_sha256=digest, expected=len(grouped),
        policy='Distinct paired heavy-atom chemical mapping classes at the lowest saved detected model count across both directions',
        hydrogen_policy='One full valid completion suffices; hydrogen permutations do not multiply choices',
        resources=dict(seconds=args.seconds, task_seconds=args.task_seconds,
                       max_mappings=args.max_mappings, workers=args.workers, memory_gib=args.memory_gib), results=completed)
    if previous is not None:
        report['previous_analysis'] = dict(sha256=sources[str(args.remaining_from.resolve())],
            runner_sha256=previous['runner_sha256'], carried_complete=len(completed),
            elapsed_seconds=previous['elapsed_seconds']+previous.get('previous_analysis',{}).get('elapsed_seconds',0))
    store = BundleStore(args.output/'bundles')
    factory = WorkPool.slurm if args.slurm else WorkPool.local
    started = time.perf_counter()
    with factory(workers=args.workers, store=store, memory_bytes=int(args.memory_gib*2**30),
                 case_seconds=args.task_seconds, shared_cache=False) as pool:
        report['ready_workers'] = len(pool.queue.ready_workers())
        carried = {r['case'] for r in completed}
        waiting = iter((case,rows) for case,rows in sorted(grouped.items()) if case not in carried)
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
                row['attempt_audit'] = asdict(pool.queue.archive_case(key,store))
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


def publish_summary(path, output):
    """Sanitized evidence and the decision-count distribution, without raw mappings."""
    import gzip
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report = json.loads(Path(path).read_bytes())
    rows = [{k:v for k,v in row.items() if k not in ('evidence','error','attempt_audit')}
            for row in report['results']]
    if len({row['case'] for row in rows}) != report['expected']:
        raise ValueError('Wait for the full analysis scope before publishing a distribution')
    exact = [row for row in rows if row['status'] == 'complete']
    pending = [row for row in rows if row['status'] != 'complete']
    counts = np.array([row['choices'] for row in exact])
    probabilities = [row['random_success'] for row in exact]
    lower_counts = sorted([*map(int,counts),*(max(1,row.get('choices_lower_bound',0)) for row in pending)])
    upper_counts = sorted([*map(int,counts),*([float('inf')]*len(pending))])
    def quantile_bounds(fraction):
        from math import ceil,isfinite
        position = ceil(fraction*len(rows))-1
        upper = upper_counts[position]
        return [lower_counts[position],int(upper) if isfinite(upper) else None]
    reference_found = sum(row.get('reference_at_minimum') is True for row in rows)
    reference_unknown = sum(row.get('reference_at_minimum') is None for row in rows)
    summary = dict(expected=report['expected'], complete_counts=len(exact), incomplete_counts=len(pending),
        reference_at_minimum=sum(row.get('reference_at_minimum') is True for row in rows),
        reference_absent_at_minimum=sum(row.get('reference_at_minimum') is False for row in rows),
        reference_unknown=sum(row.get('reference_at_minimum') is None for row in rows),
        detected_minimum_proven=sum(row.get('saved_minimum_proven',False) for row in rows),
        exact_count_distribution=dict(Counter(map(int,counts))),
        single_choice=int(np.count_nonzero(counts == 1)),
        exact_subset_mean=float(counts.mean()) if len(counts) else None,
        exact_subset_median=float(np.median(counts)) if len(counts) else None,
        exact_subset_p95=float(np.percentile(counts,95)) if len(counts) else None,
        exact_subset_max=int(max(counts)) if len(counts) else None,
        exact_subset_random_accuracy=float(np.mean(probabilities)) if probabilities else None,
        full_scope_median_bounds=quantile_bounds(.5),full_scope_p95_bounds=quantile_bounds(.95),
        full_scope_reference_at_minimum_bounds=[reference_found/len(rows),
                                               (reference_found+reference_unknown)/len(rows)],
        full_scope_random_accuracy_bounds=[sum(probabilities)/len(rows),
            (sum(probabilities)+sum(1/max(1,row.get('choices_lower_bound',0)) for row in pending))/len(rows)],
        analysis_wall_seconds=report['elapsed_seconds']+report.get('previous_analysis',{}).get('elapsed_seconds',0),
        ready_workers=report['ready_workers'], source_unchanged=report['source_unchanged'],
        policy=report['policy'], hydrogen_policy=report['hydrogen_policy'], resources=report['resources'],
        runner_sha256=report['runner_sha256'], report_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest())
    output = Path(output); output.mkdir(parents=True,exist_ok=True)
    _publish(output/'summary.json',summary)
    with gzip.open(output/'reactions.json.gz','wt') as stream:
        json.dump(sorted(rows,key=lambda row:row['case']),stream)
    edges = [1,2,5,10,20,50,100,float('inf')]
    labels = ['1','2','3–5','6–10','11–20','21–50','51–100','>100']
    groups = np.zeros((2,len(edges)),dtype=int)
    for row in exact:
        i = next(i for i,edge in enumerate(edges) if row['choices'] <= edge)
        groups[int(not row['reference_at_minimum']),i] += 1
    used = max(np.flatnonzero(groups.sum(axis=0)),default=0)+1
    labels,groups = labels[:used],groups[:,:used]
    fig,ax = plt.subplots(figsize=(8,4.5),layout='constrained')
    x=np.arange(len(labels))
    ax.bar(x,groups[0],color='#287b8e',label='Reference among these options')
    ax.bar(x,groups[1],bottom=groups[0],color='#d18a59',label='Reference absent')
    for i,total in enumerate(groups.sum(axis=0)):
        if total:
            ax.annotate(f'{total:,}',(i,total),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
    ax.set(xticks=x,xticklabels=labels,xlabel='Distinct mapping options at the detected minimum',
           ylabel='Golden reactions (log scale)',yscale='log',title='Choices remaining under minimum-event ranking')
    ax.set_ylim(.8,max(2,int(groups.sum(axis=0).max())*2))
    ax.spines[['top','right']].set_visible(False)
    ax.legend(frameon=False,fontsize=9,loc='upper right')
    fig.supxlabel(f'Exact counts: {len(exact):,}/{len(rows):,} reactions; incomplete counts: {len(pending):,}. '
                  'Equivalent heavy mappings and hydrogen shuffles count once.',fontsize=9)
    for suffix in ('png','pdf','svg'):
        fig.savefig(output/('minimum-choices.'+suffix),dpi=180)
    plt.close(fig)
    return summary


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
    parser.add_argument('--remaining-from',type=Path,help='Preserve exact counts and analyze only incomplete reactions')
    parser.add_argument('--plot', action='store_true',help='Summarize saved results; campaign argument is results.json')
    args = parser.parse_args()
    if args.plot:
        print(json.dumps(publish_summary(args.campaign,args.output),indent=2))
    else:
        run(args)
