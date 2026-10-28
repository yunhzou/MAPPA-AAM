"""Prove Golden event minima from saved families with its original score policy.

This stage does not enumerate output event classes. It establishes per-seed
minima and valid decoding bounds, including partial mappings and both directions.
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import tarfile
import tempfile
import time
import zipfile

from mappa.execution import Task
from .saved_analysis import prepare_cut, read_member


def plan_minima(archive, task_id, *, seconds=3600, memory_bytes=4*2**30):
    archive = Path(archive).resolve()
    with tarfile.open(archive) as tar:
        manifest = json.load(tar.extractfile('archive-manifest.json'))['files']
        def member(name):
            item = tar.getmember(name)
            assert item.isfile() and item.size == manifest[name]['bytes'], name
            return dict(archive=str(archive), name=name, offset=item.offset_data, **manifest[name])
        source = {name: member(name) for name in ('input.json', 'reference.json', 'task.json')}
        search = json.loads(read_member(member('search.json')))
        tasks = []
        for cut in range(search['metrics']['cut_count']):
            raw = f'cuts/cut_{cut:05d}.raw.pkl.gz'
            final = f'finalized/{cut:05d}.pkl.gz'
            graph = final if final in manifest else raw
            payload = dict(source=source, graph=member(graph), raw=member(raw),
                finalized=final in manifest, cut=cut, task_id=task_id,
                family_seconds=30., cut_seconds=3000.)
            tasks.append(Task(f'{task_id:05d}/minimum/{cut:05d}',
                'rebenchmark.minima:minimum_cut', payload, case=f'{task_id:05d}',
                case_seconds=seconds, memory_bytes=memory_bytes))
        return tasks


def family_key(path):
    """Equality covers every input consumed by the existing path compiler."""
    from mappa.search_graph import frozen_value
    placements = []
    for edge in path.transitions:
        fragment = path.graph.fragment_placement(edge)
        if fragment is not None:
            locked = path.graph.states[path.graph.transitions[edge].source].mapping
            placements.append((locked, frozen_value(asdict(fragment))))
    return (tuple(sorted(path.mapping.items())), path.context.cuts,
            path.context.graph_floor, path.context.iso_tolerance, tuple(placements))


def selected_families(plan, graph):
    """Use Golden's best heavy/all-atom coverage tier independently per seed."""
    tiers = {}
    for terminal in graph.terminals:
        state = graph.states[terminal]
        mapping = dict(state.mapping)
        tier = (-sum(plan.problem.reactant.elements[r] != 'H' for r in mapping), -len(mapping))
        tiers[state.context] = min(tier, tiers.get(state.context, tier))
    families = {}
    for terminal in graph.terminals:
        state = graph.states[terminal]
        mapping = dict(state.mapping)
        tier = (-sum(plan.problem.reactant.elements[r] != 'H' for r in mapping), -len(mapping))
        if tier != tiers[state.context]:
            continue
        for path in graph.paths(terminal):
            key = family_key(path)
            if key not in families:
                families[key] = [path, []]
            families[key][1].append(dict(seed=state.context+1, terminal=terminal,
                                        transitions=path.transitions))
    return tiers, list(families.values())


def seed_bounds(tiers, families, results):
    rows = {}
    for seed, tier in tiers.items():
        selected = [i for i, (_, origins) in enumerate(families)
                    if any(p['seed'] == seed+1 for p in origins)]
        completed = [results[i] for i in selected if i in results]
        upper = min((r['upper_bound'] for r in completed), default=None)
        lower = min((results[i]['lower_bound'] if i in results else 0 for i in selected), default=0)
        # A tight global bound is sufficient even if a worse family is unfinished.
        rows[str(seed+1)] = dict(coverage_tier=tier, lower_bound=lower, upper_bound=upper,
            minimum_proven=upper is not None and lower == upper,
            families=len(selected), processed_families=len(completed),
            unfinished_families=[i for i in selected if i not in results or not results[i]['optimal']])
    return rows


def minimum_cut(payload, workspace):
    from mappa.family_scoring import minimize_events
    plan, reference, graph, graph_digest, finalized = prepare_cut(payload, workspace)
    started = time.monotonic()
    deadline = started+payload['cut_seconds']
    tiers, families = selected_families(plan, graph)
    results = {}
    with tempfile.TemporaryDirectory(dir=workspace.scratch) as td:
        evidence = Path(td)/'minimum-evidence.zip'
        with zipfile.ZipFile(evidence, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for number, (path, origins) in enumerate(families):
                left = deadline-time.monotonic()
                if left <= 0:
                    break
                queries = []
                def record_query(solver, status, bound):
                    name = f'{number:06d}/{len(queries):04d}.smt2'
                    text = solver.to_smt2()
                    archive.writestr(name, text)
                    queries.append(dict(name=name, status=status, bound=bound,
                                        sha256=hashlib.sha256(text.encode()).hexdigest()))
                result = minimize_events(path, plan.problem, reverse=plan.reversed,
                    seconds=min(payload['family_seconds'], left), record_query=record_query)
                results[number] = result
                archive.writestr(f'{number:06d}/result.json', json.dumps(dict(result=result,
                    provenance=origins, queries=queries)))
        report = dict(task_id=payload['task_id'], cut=payload['cut'],
            seeds=seed_bounds(tiers, families, results), family_count=len(families),
            processed_families=len(results), wall_seconds=time.monotonic()-started,
            policy=dict(name='original_golden', bond_floor=.2, event_tolerance=.5,
                        coverage='best heavy/all-atom coverage tier', reverse=plan.reversed),
            source=payload['source'], raw=payload['raw'], graph=payload['graph'],
            finalized_sha256=graph_digest, finalized=finalized,
            evidence=workspace.writer.put_file(evidence), search_rerun=False,
            scope='Proven minima/bounds only; not completed event-class enumeration')
        return report
