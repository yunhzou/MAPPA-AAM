"""Continue missing cut/seed checks from verified saved search archives.

The shared executor handles CPU placement and the directional wall budget.
Archive byte ranges avoid copying old solver logs or extracting shared trees.
"""
import hashlib
import json
from pathlib import Path
import tarfile
import tempfile

from mappa.mappa_core.core.scheduler import Task
from .evidence import read


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
    from .inputs import load_plan
    from mappa.mappa_core.utils.io.checkpoints import read_graph_checkpoint, read_raw_cut, write_graph_checkpoint
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
                from mappa.mappa_core.core.matching.graphs import build_graph
                from mappa.mappa_core.core.matching.symmetry.workspace import ConditionedSymmetryWorkspace
                from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry
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
