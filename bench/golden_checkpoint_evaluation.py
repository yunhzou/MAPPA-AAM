"""Reference-family verification over independently persisted search cuts."""
import gc
import hashlib
import json
from pathlib import Path
import time

from mappa.aam import checkpoint_manifest
from mappa.alignment.sweep import cut_sweep_items
from mappa.artifacts import raw_cut_paths, read_raw_cut
from mappa.conditioned_symmetry import ConditionedSymmetryWorkspace
from mappa.domain import AAMResult, AAMSearchMetrics
from mappa.frag import build_graph
from mappa.search_symmetry import finalize_graph_symmetry
from golden_evaluation import evaluate_planned


def aggregate_cut_verdicts(expected_count, rows):
    """An exact negative requires every scheduled cut, not just every saved cut."""
    if expected_count < 1:
        raise ValueError('A cut schedule must include its uncut run')
    indices = [row['cut_index'] for row in rows]
    if len(indices) != len(set(indices)) or any(i not in range(expected_count) for i in indices):
        raise ValueError('Duplicate or unexpected cut verdict')
    outcomes = [row['reference_recovery'] for row in rows]
    if any(s not in ('recovered', 'not_recovered', 'unknown') for s in outcomes):
        raise ValueError('Invalid cut verdict')
    if 'recovered' in outcomes:
        return 'recovered'
    if len(indices) == expected_count and all(s == 'not_recovered' for s in outcomes):
        return 'not_recovered'
    return 'unknown'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def evaluate_checkpoints(directory, plan, features, reference, *, seconds=275,
                         cut_seconds=20, query_timeout_ms=5000, progress=None):
    """Check one finalized cut at a time; returned witness/ranking is cut-local.

    Independent cut graphs form a disjoint union. Definitive negative results
    are sound only for a complete scheduled union with matching input identity.
    The original search can have failed while combining already-complete cuts.
    """
    directory = Path(directory)
    config, problem = plan.config, plan.problem
    expected_manifest = checkpoint_manifest(problem, config)
    manifest = json.loads((directory / 'manifest.json').read_text())
    if manifest.get('execution') == 'shared_policies':
        expected_manifest['execution'] = 'shared_policies'
    if manifest != expected_manifest:
        raise ValueError('Checkpoint input or configuration differs from verification')
    cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor) if config.sweep_cuts else [()]
    files = raw_cut_paths(directory)
    indices = [int(path.name.split('_')[1].split('.')[0]) for path in files]
    if any(i not in range(len(cuts)) for i in indices):
        raise ValueError('Unexpected cut checkpoint index')
    started = time.perf_counter()
    cpu = time.process_time()
    deadline = started + seconds
    target = build_graph(problem.product.elements, problem.product.wbo, bond_cut=config.graph_floor)
    rows = []
    witness = None
    for index, path in zip(indices, files):
        if time.perf_counter() >= deadline:
            break
        before = time.perf_counter()
        graph = read_raw_cut(path)
        if not graph.contexts or any(tuple(c.cuts) != tuple(cuts[index]) for c in graph.contexts):
            raise ValueError('Saved graph has the wrong cut context')
        # A fresh workspace prevents conditioned-transition caches accumulating
        # across cuts. Full matching histories remain in the raw checkpoints.
        workspace = ConditionedSymmetryWorkspace(target, config.iso_tolerance)
        graph, _ = finalize_graph_symmetry(graph, target,
            iso_tolerance=config.iso_tolerance, workspace=workspace)
        del workspace
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            del graph
            break
        result = evaluate_planned(AAMResult(problem, config, graph, AAMSearchMetrics.from_record({}, 0)),
            plan, features, reference, seconds=min(cut_seconds, remaining),
            query_timeout_ms=query_timeout_ms)
        rows.append(dict(cut_index=index, cut=str(path), sha256=sha(path),
            reference_recovery=result['reference_recovery'],
            wall_seconds=time.perf_counter() - before))
        del graph
        if result['reference_recovery'] == 'recovered':
            witness = dict(cut=str(path), cut_index=index, evaluation=result)
        else:
            del result
        gc.collect()
        if progress is not None:
            progress(dict(reference_recovery=aggregate_cut_verdicts(len(cuts), rows),
                          expected_cuts=len(cuts), available_cuts=len(files), checked=list(rows)))
        if witness is not None:
            break
    verdict = aggregate_cut_verdicts(len(cuts), rows)
    return dict(reference_recovery=verdict, witness=witness,
        evaluation_scope='Reference-family membership in the saved cut union; witness and ranking, if present, are cut-local.',
        checkpoint_verification=dict(expected_cuts=len(cuts), available_cuts=len(files), checked=rows,
            all_expected_cuts_checked={r['cut_index'] for r in rows} == set(range(len(cuts))),
            manifest_sha256=sha(directory / 'manifest.json'), verifier_sha256=sha(__file__)),
        cpu_seconds=time.process_time() - cpu, wall_seconds=time.perf_counter() - started)
