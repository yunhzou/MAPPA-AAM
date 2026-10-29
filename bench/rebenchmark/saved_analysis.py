"""Continue missing cut/seed checks from verified saved search archives.

The shared executor handles CPU placement and the directional wall budget.
Archive byte ranges avoid copying old solver logs or extracting shared trees.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import secrets
import socket
import tarfile
import tempfile
import time

from mappa.execution import Task, TaskQueue, coordinator, connect, run, read_blob
from .evidence import read, save


def read_member(member):
    with open(member['archive'], 'rb') as stream:
        stream.seek(member['offset'])
        data = stream.read(member['bytes'])
    if len(data) != member['bytes'] or hashlib.sha256(data).hexdigest() != member['sha256']:
        raise ValueError('Saved archive member checksum differs: '+member['name'])
    return data


def plan_missing(archive, task_id, *, seconds=3600, memory_bytes=4*2**30):
    """Completed evaluations are immutable; only absent records become tasks."""
    archive = Path(archive).resolve()
    with tarfile.open(archive) as tar:
        manifest = json.load(tar.extractfile('archive-manifest.json'))['files']
        def member(name):
            entry = tar.getmember(name)
            if not entry.isfile() or entry.size != manifest[name]['bytes']:
                raise ValueError('Saved archive member size differs: '+name)
            return dict(archive=str(archive), name=name, offset=entry.offset_data,
                        **manifest[name])
        source = {name: member(name) for name in ('input.json', 'reference.json', 'task.json')}
        task = json.loads(read_member(source['task.json']))
        search = json.loads(read_member(member('search.json')))
        count = search['metrics']['cut_count']
        tasks = []
        for cut in range(count):
            final_name = f'finalized/{cut:05d}.pkl.gz'
            raw_name = f'cuts/cut_{cut:05d}.raw.pkl.gz'
            for seed in range(task['config']['seed_count']):
                record = f'analysis/{cut:05d}/{seed:02d}/evaluation.json'
                if record in manifest:
                    continue
                graph_name = final_name if final_name in manifest else raw_name
                payload = dict(source=source, graph=member(graph_name), raw=member(raw_name),
                    finalized=final_name in manifest, cut=cut, seed=seed, task_id=task_id)
                tasks.append(Task(f'{task_id:05d}/cut/{cut:05d}/seed/{seed:02d}',
                    'rebenchmark.saved_analysis:evaluate', payload, case=f'{task_id:05d}',
                    case_seconds=seconds, memory_bytes=memory_bytes))
        return tasks


def prepare_cut(payload, workspace):
    from .campaign import load_plan
    from mappa.artifacts import read_graph_checkpoint, read_raw_cut, write_graph_checkpoint
    def prepare():
        with tempfile.TemporaryDirectory(dir=workspace.scratch) as td:
            folder = Path(td)
            for name, member in payload['source'].items():
                (folder/name).write_bytes(read_member(member))
            plan = load_plan(folder)
            graph_path = folder/Path(payload['graph']['name']).name
            graph_path.write_bytes(read_member(payload['graph']))
            if payload['finalized']:
                graph = read_graph_checkpoint(graph_path)
                digest = payload['graph']['sha256']
                finalized = None
            else:
                from mappa.frag import build_graph
                from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
                from mappa.mappa_core.core.search.symmetry import finalize_graph_symmetry
                raw = read_raw_cut(graph_path)
                target = build_graph(plan.problem.product.elements, plan.problem.product.wbo,
                                     bond_cut=plan.config.graph_floor)
                graph, _ = finalize_graph_symmetry(raw, target, iso_tolerance=plan.config.iso_tolerance,
                    workspace=ConditionedSymmetryWorkspace(target, plan.config.iso_tolerance))
                path = folder/'finalized.pkl.gz'
                write_graph_checkpoint(graph, path)
                data = path.read_bytes()
                digest = hashlib.sha256(data).hexdigest()
                finalized = workspace.writer.put_file(path)
            return plan, read(folder/'reference.json'), graph, digest, finalized
    key = ('saved-cut', payload['graph']['archive'], payload['graph']['sha256'],
           payload['source']['task.json']['sha256'])
    return workspace.prepared(key, prepare)


def evaluate(payload, workspace):
    from .analysis import analyze_seed
    plan, reference, graph, digest, finalized = prepare_cut(payload, workspace)
    with tempfile.TemporaryDirectory(dir=workspace.scratch) as td:
        folder = Path(td)/'evidence'
        summary = analyze_seed(plan, reference, graph, payload['seed'], folder, payload['cut'],
            raw_sha256=payload['raw']['sha256'], finalized_sha256=digest)
        archive = Path(td)/'evidence.tar.gz'
        with tarfile.open(archive, 'w:gz') as tar:
            tar.add(folder, arcname='evidence')
        return dict(task_id=payload['task_id'], cut=payload['cut'], seed=payload['seed']+1,
            evaluation=read(folder/'evaluation.json'), classes=read(folder/'classes.json'),
            summary=summary, evidence=workspace.writer.put_file(archive),
            source=payload['source'], raw=payload['raw'], graph=payload['graph'],
            finalized=finalized, search_rerun=False)


def execute(args):
    output = args.output
    specification = read(output/'plan.json')
    memory = specification['memory_bytes']
    identity = specification['identity']
    if args.mode == 'worker':
        endpoint = read(output/'endpoint.json')
        queue = connect((endpoint['host'], endpoint['port']), (output/'worker.key').read_bytes())
        run(queue, workers=args.workers, identity=identity, scratch=os.environ['TMPDIR'],
            bundles=output/'bundles', memory_bytes=memory)
        return
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        queue = TaskQueue(Path(scratch)/'queue.sqlite', identity)
        tasks = []
        for row in specification['archives']:
            tasks.extend(plan_missing(row['path'], row['task'],
                seconds=specification['case_seconds'], memory_bytes=memory))
        if specification.get('selected_keys') is not None:
            selected = set(specification['selected_keys'])
            tasks = [task for task in tasks if task.key in selected]
            assert {task.key for task in tasks} == selected
        queue.add(tasks)
        save(output/'tasks.json', [asdict(task) for task in tasks])
        secret = secrets.token_bytes(32)
        with (output/'worker.key').open('xb') as stream:
            os.chmod(output/'worker.key', 0o600)
            stream.write(secret)
        with coordinator(queue, (socket.gethostname(), 0), secret) as address:
            save(output/'endpoint.json', dict(host=address[0], port=address[1], tasks=len(tasks)))
            # This process owns the queue; separate allocations execute tasks.
            last_backup = 0.
            while any(n for state, n in queue.snapshot().items() if state in ('running', 'pending')):
                if time.monotonic()-last_backup >= 60:
                    queue.backup(output/'queue.sqlite')
                    last_backup = time.monotonic()
                save(output/'progress.json', dict(states=queue.snapshot(), utc=time.time()))
                time.sleep(10)
            queue.backup(output/'queue.sqlite')
            rows = queue.results()
            merged = {}
            for key, row in rows.items():
                record = dict(state=row['state'], resource=row['result'])
                if row['state'] == 'complete':
                    value = read_blob(row['result']['blob'], scratch)
                    record['evaluation'] = value['evaluation']
                    record['classes'] = value['classes']
                    record['evidence'] = asdict(value['evidence'])
                    record['source'] = value['source']
                    record['raw'] = value['raw']
                    record['graph'] = value['graph']
                    record['finalized'] = asdict(value['finalized']) if value['finalized'] else None
                if record['resource'] and 'blob' in record['resource']:
                    record['resource'] = dict(record['resource'], blob=asdict(record['resource']['blob']))
                merged[key] = record
            save(output/'results.json', merged)
            save(output/'progress.json', dict(states=queue.snapshot(), utc=time.time(), finished=True))
            # Let remote workers receive the terminal queue state before closing.
            time.sleep(5)
        queue.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('coordinator', 'worker'))
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    execute(parser.parse_args())
