"""Two-node execution validation, not a new Golden benchmark campaign.

Run one rank per node under Slurm. All scientific controls are computed from
unchanged public search/evaluation operations; the pilot compares execution only.
"""
import argparse
from collections import Counter
from dataclasses import replace
import json
import os
from pathlib import Path
import secrets
import socket
import tempfile
import time

from mappa import AAMProblem, AAMSearchConfig, AAMSearchPlan
from mappa.aam import SearchSession
from mappa.alignment.sweep import cut_sweep_items
from mappa.execution import TaskQueue, BundleWriter, coordinator, connect, run, read_blob
from mappa.execution.search import iter_cuts
from mappa.search_symmetry import finalize_graph_symmetry
from golden_evaluation import prepare
from .shared import plan_direction
from .analysis import analyze_seed
from .evidence import save, read


def wait_file(path, seconds=240):
    deadline = time.monotonic()+seconds
    while not path.exists():
        if time.monotonic() > deadline:
            raise TimeoutError(str(path))
        time.sleep(.2)
    return read(path)


def controls():
    reactions = [
        '[CH3:1][OH:2]>>[CH3:1][OH:2].[Cl-]',
        '[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][OH:3]',
    ]
    config = AAMSearchConfig(seed_count=10, branch_limit=100, iso_tolerance=1)
    for index, reaction in enumerate(reactions):
        original, features, mapping = prepare(reaction)
        for reverse in (False, True):
            directed = AAMProblem(original.product, original.reactant, original.name) if reverse else original
            plan = AAMSearchPlan(original, directed, config, reverse)
            yield f'control-{index}-{plan.direction}', plan, dict(features=features, mapping=sorted(mapping.items()))


def validate(queue, cases, scratch):
    rows, checks = queue.results(), Counter()
    cpu, wall = time.process_time(), time.perf_counter()
    for case, plan, reference, cut_keys, evaluation_keys in cases:
        session = SearchSession(plan.problem, plan.config)
        for number, (cut, observed) in enumerate(zip(cut_sweep_items(plan.problem.reactant.wbo, plan.config.cut_floor),
                                                    iter_cuts(queue, cut_keys, scratch))):
            graph, _ = session.search_cut(cut)
            graph, _ = finalize_graph_symmetry(graph, session.context[2], iso_tolerance=plan.config.iso_tolerance)
            assert graph == observed['graph'], (case, number)
            checks['equal_cut_graphs'] += 1
            for seed in range(len(graph.contexts)):
                key = evaluation_keys[number*plan.config.seed_count+seed]
                actual = read_blob(rows[key]['result']['blob'], scratch)
                with tempfile.TemporaryDirectory(dir=scratch) as temporary:
                    folder = Path(temporary)
                    expected = analyze_seed(plan, reference, graph, seed, folder, number)
                    assert expected == actual['summary'], (case, number, seed)
                    assert read(folder/'classes.json') == actual['classes']
                    verdict = read(folder/'evaluation.json')
                    assert verdict['reference_recovery'] == actual['evaluation']['reference_recovery']
                checks['equal_seed_evaluations'] += 1
    return dict(checks=checks, direct_comparison_cpu_seconds=time.process_time()-cpu,
                direct_comparison_wall_seconds=time.perf_counter()-wall)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--identity', required=True)
    args = parser.parse_args()
    rank = int(os.environ['SLURM_PROCID'])
    assert int(os.environ['SLURM_NTASKS']) == 2
    output = args.output
    output.mkdir(exist_ok=True, parents=True)
    for batch in (1, 2, 5, 10):
        folder = output/f'batch-{batch}'
        if rank == 0:
            folder.mkdir(exist_ok=True)
            with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
                queue = TaskQueue(Path(scratch)/'queue.db', args.identity)
                writer = BundleWriter(folder/'inputs', scratch)
                cases = []
                for name, plan, reference in controls():
                    tasks, cuts, evaluations = plan_direction(name, plan, reference, writer, seed_batch=batch,
                                                              case_seconds=600, memory_bytes=1024*2**20)
                    queue.add(tasks)
                    cases.append((name, plan, reference, cuts, evaluations))
                secret = secrets.token_bytes(32)
                auth = folder/'worker.key'
                with auth.open('xb') as stream:
                    os.chmod(auth, 0o600)
                    stream.write(secret)
                with coordinator(queue, (socket.gethostname(), 0), secret) as address:
                    save(folder/'endpoint.json', dict(host=address[0], port=address[1]))
                    start = time.perf_counter()
                    run(queue, workers=3, identity=args.identity, scratch=scratch,
                        bundles=folder/'results', memory_bytes=1024*2**20)
                    wait_file(folder/'rank-1.json', seconds=120)
                    elapsed = time.perf_counter()-start
                    assert set(queue.snapshot()) == {'complete'}, queue.results()
                    rows = queue.results()
                    workers = Counter(row['result']['worker'].split(':')[0] for row in rows.values())
                    assert len(workers) == 2, workers
                    queue.backup(folder/'queue-backup.sqlite')
                verification = validate(queue, cases, scratch)
                measures = dict(batch=batch, wall_seconds=elapsed, tasks=len(rows), tasks_by_host=workers,
                    compute_cpu_seconds=sum(r['result']['cpu_seconds'] for r in rows.values()),
                    compute_seconds=sum(r['result']['compute_seconds'] for r in rows.values()),
                    publish_seconds=sum(r['result']['publish_seconds'] for r in rows.values()),
                    monitoring_seconds=sum(r['result']['monitoring_seconds'] for r in rows.values()),
                    max_worker_rss_bytes=max(r['result']['peak_tree_rss_bytes'] for r in rows.values()),
                    shared_files=sum(1 for p in folder.rglob('*') if p.is_file()),
                    **verification)
                save(folder/'validation.json', measures)
                print(json.dumps(measures), flush=True)
                queue.close()
                auth.unlink()
        else:
            endpoint = wait_file(folder/'endpoint.json')
            key = (folder/'worker.key').read_bytes()
            queue = connect((endpoint['host'], endpoint['port']), key)
            state = run(queue, workers=3, identity=args.identity, scratch=os.environ['TMPDIR'],
                        bundles=folder/'results', memory_bytes=1024*2**20)
            save(folder/'rank-1.json', dict(host=socket.gethostname(), state=state))
    if rank == 0:
        save(output/'validation.json', dict(complete=True, scope='Two small reactions, both directions; execution controls only',
             batches=[read(output/f'batch-{batch}'/'validation.json') for batch in (1, 2, 5, 10)]))


if __name__ == '__main__':
    main()
