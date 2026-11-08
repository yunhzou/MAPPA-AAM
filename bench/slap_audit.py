"""Independent SLAP result audit and its standard, uncut binary baseline."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import random
import sqlite3
import tempfile
import time

from compare_slap import hardware, slap_source
from golden import _measure, _publish
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob


def paired_graph(features, mapping):
    """Direct molecular edges plus cross-endpoint mapping edges, without nauty."""
    import networkx as nx
    graph = nx.Graph()
    for side, feature in enumerate(features):
        for atom, color in zip(feature['heavy'], feature['colors']):
            graph.add_node((side, atom), color=(side, *color))
        for a, b, color in feature['bonds']:
            graph.add_edge((side, feature['heavy'][a]), (side, feature['heavy'][b]),
                           color=('bond', *color))
    sources, targets = set(), set()
    for a, b in mapping:
        if a in sources or b in targets or (0, a) not in graph or (1, b) not in graph:
            raise ValueError('Mapping is not an injective heavy-atom correspondence')
        sources.add(a); targets.add(b)
        graph.add_edge((0, a), (1, b), color=('mapping',))
    return graph


def equivalent(left, right):
    """Independent VF2 oracle; return one joint endpoint correspondence."""
    from networkx.algorithms.isomorphism import GraphMatcher
    matcher = GraphMatcher(left, right,
        node_match=lambda a, b: a['color'] == b['color'],
        edge_match=lambda a, b: a['color'] == b['color'])
    return next(matcher.isomorphisms_iter(), None)


def label_mapping(labels, elements):
    """Interpret native label groups, validating balancing atoms as well."""
    from rdkit import Chem
    table = Chem.GetPeriodicTable()
    numbers = [[table.GetAtomicNumber(str(e)) for e in side] for side in elements]
    difference = Counter(numbers[1]); difference.subtract(numbers[0])
    for element, count in difference.items():
        numbers[0 if count > 0 else 1].extend([element] * abs(count))
    groups = []
    for side in range(2):
        if len(labels[side]) != len(numbers[side]):
            raise ValueError('Native labels do not cover the balanced input')
        grouped = defaultdict(list)
        for atom, label in enumerate(labels[side]):
            grouped[int(label)].append(atom)
        groups.append(grouped)
    if set(groups[0]) != set(groups[1]):
        raise ValueError('Native endpoint labels differ')
    pairs = []
    for label, source in groups[0].items():
        target = groups[1][label]
        kinds = {numbers[0][i] for i in source} | {numbers[1][i] for i in target}
        if len(source) != len(target) or len(kinds) != 1:
            raise ValueError('Native label group violates element or cardinality')
        if kinds == {1}:
            continue
        if len(source) != 1:
            raise ValueError('Heavy label is not a unique native assignment')
        a, b = source[0], target[0]
        if a < len(elements[0]) and b < len(elements[1]):
            pairs.append((a, b))
    return tuple(sorted(pairs))


def audit_saved(payload, workspace):
    from slapmapper._idxmapstr import lgp2idx_map_str, parse_index_mapping_string
    from types import SimpleNamespace
    raw = read_blob(Blob(**payload['input']), workspace.scratch)
    reference = read_blob(Blob(**payload['reference']), workspace.scratch)
    matched = workspace.load(Blob(**payload['matching']))
    saved = workspace.load(Blob(**matched['artifact']))
    verified = workspace.load(Blob(**payload['verification']))['verification']
    features = reference['features'][::-1] if payload['reversed'] else reference['features']
    truth = [(b, a) for a, b in reference['mapping']] if payload['reversed'] else reference['mapping']
    elements = [raw[side]['elements'] for side in ['reactant', 'product']]
    mappings, uncut, calls = set(), set(), set()
    export_checks = 0
    for record in saved['records']:
        mapping = label_mapping(record['labels'], elements)
        if mapping != tuple(map(tuple, record['mapping'])):
            raise ValueError('Saved pairs differ from native labels')
        mappings.add(mapping)
        calls.add((record['mode'], tuple(map(tuple, record['cut']))))
        if not record['cut']:
            # Independently check the public index-string export on uncut outputs.
            graphs = []
            for labels in record['labels']:
                groups = defaultdict(list)
                for atom, label in enumerate(labels):
                    groups[label].append(atom)
                graphs.append(SimpleNamespace(label2idxs=groups))
            groups = parse_index_mapping_string(lgp2idx_map_str(graphs))
            exported = tuple(sorted((a[0], b[0]) for a, b in groups
                if len(a) == len(b) == 1 and a[0] < len(elements[0])
                and b[0] < len(elements[1]) and elements[0][a[0]] != 'H'))
            if exported != mapping:
                raise ValueError('Public SLAP export differs from saved heavy mapping')
            export_checks += 1
            if record['mode'] == 'binary':
                uncut.add(mapping)
    if mappings != set(tuple(map(tuple, m)) for m in saved['mappings']):
        raise ValueError('Saved mapping union differs from native records')
    if len(calls) != matched['search']['calls']:
        raise ValueError('A completed mode/cut call is missing')
    expected = paired_graph(features, truth)
    candidates = [tuple(map(tuple, verified['witness']))] if verified['status'] == 'recovered' else sorted(mappings)
    if not set(candidates) <= mappings:
        raise ValueError('Stored verification witness is not a returned mapping')
    checked = 0
    for mapping in candidates:
        checked += 1
        if equivalent(expected, paired_graph(features, mapping)) is not None:
            status = 'recovered'; break
    else:
        status = 'not_recovered'
    if status != verified['status']:
        raise ValueError('Independent VF2 verdict disagrees with saved verification')
    binary_status = 'recovered' if any(equivalent(expected, paired_graph(features, m))
        is not None for m in sorted(uncut)) else 'not_recovered'
    return dict(status=status, native_records=len(saved['records']), mappings=len(mappings),
        oracle_checks=checked, export_checks=export_checks, calls=len(calls),
        binary_uncut_status=binary_status, cpu_model=hardware())


def native_call(payload, workspace):
    """Time default binary SLAP on the identical prepared chemical endpoints."""
    import numpy as np
    from slapmapper.aam import SlapAAM
    from compare_slap import slap_graphs
    from mappa.mappa_core.utils.io.records import problem_from_record
    raw = read_blob(Blob(**payload['input']), workspace.scratch)
    reference = read_blob(Blob(**payload['reference']), workspace.scratch)
    problem = problem_from_record(raw)
    features = reference['features'][::-1] if payload['reversed'] else reference['features']
    elements = [raw[side]['elements'] for side in ['reactant', 'product']]
    random.seed(20260910); np.random.seed(20260910)
    metrics = {}
    with _measure(metrics, 'standard_search'):
        graphs = slap_graphs(problem)
        targets = [i for i, z in enumerate(graphs[0].labels) if z > 1]
        mapper = SlapAAM()
        mapper.get_maps(graphs, break_sym_targets=targets)
        mappings = {label_mapping([g.labels for g in r['lgp']], elements) for r in mapper.results}
    truth = [(b, a) for a, b in reference['mapping']] if payload['reversed'] else reference['mapping']
    expected = paired_graph(features, truth)
    hit = correspondence = None
    with _measure(metrics, 'verification'):
        for mapping in sorted(mappings):
            correspondence = equivalent(expected, paired_graph(features, mapping))
            if correspondence is not None:
                hit = mapping
                break
    return dict(status='recovered' if hit is not None else 'not_recovered',
                mappings=len(mappings), candidates=sorted(mappings), witness=hit,
                endpoint_correspondence=list(correspondence.items()) if correspondence is not None else None,
                metrics=metrics, cpu_model=hardware())


def run(args):
    from mappa.execution import BundleStore, WorkPool
    from mappa.mappa_core.core.scheduler import Task
    from enzymemap import Ledger, require_drained
    source = args.campaign.resolve()
    require_drained(source)
    plan = json.loads((source/'plan.json').read_text())
    if hardware() != plan['cpu_model']:
        raise ValueError('Audit must use the paired benchmark CPU model')
    if slap_source() != plan['slap_hashes']:
        raise ValueError('Pinned SLAP source changed')
    for variable in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
        if os.environ.get(variable) != '1':
            raise ValueError('Use one numerical-library thread')
    with sqlite3.connect(f'file:{source}/ledger.sqlite?mode=ro', uri=True) as db:
        records = {(n,s):(state,ref) for n,s,state,ref in db.execute('SELECT case_name,stage,state,reference FROM results')}
        source_rows = db.execute('SELECT * FROM results ORDER BY case_name,stage').fetchall()
    source_digest = hashlib.sha256(json.dumps(source_rows,separators=(',',':')).encode()).hexdigest()
    tasks=[]
    for row in plan['inputs']:
        name=f"slap/{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
        for stage in ['native','audit']:
            payload=dict(row, name=name)
            if stage=='audit':
                if records[name,'matching'][0]!='complete':
                    continue
                payload.update(matching=json.loads(records[name,'matching'][1]),
                               verification=json.loads(records[name,'verification'][1]))
            tasks.append((name,stage,payload))
    args.output.mkdir(parents=True,exist_ok=False)
    _publish(args.output/'plan.json',dict(source_plan_sha256=hashlib.sha256((source/'plan.json').read_bytes()).hexdigest(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), tasks=len(tasks),
        cpu_model=hardware(),slap_hashes=slap_source(),memory_gib=4,seconds=600,
        source_rows_sha256=source_digest))
    with tempfile.TemporaryDirectory() as scratch:
        ledger=Ledger(args.output,scratch,resume=False)
        store=BundleStore(args.output/'bundles',scratch=scratch)
        with WorkPool.slurm(workers=args.workers,store=store,memory_bytes=4*2**30,
                            case_seconds=600,shared_cache=False) as pool:
            pending=iter(tasks);active={};exhausted=False;next_save=0.; stopped_on_error=None
            while active or not exhausted:
                while len(active)<2*args.workers and not exhausted:
                    try:name,stage,payload=next(pending)
                    except StopIteration:exhausted=True;break
                    key=name+'/'+stage
                    operation='native_call' if stage=='native' else 'audit_saved'
                    pool.queue.add([Task(key,'slap_audit:'+operation,payload,case=key,case_seconds=600)])
                    active[key]=(name,stage)
                for key,(name,stage) in list(active.items()):
                    state=pool.queue.results([key])[key]
                    if state['state'] in ['pending','running']:continue
                    value=state['result'] or {};complete=state['state']=='complete'
                    ledger.put(name,stage,'complete' if complete else 'error',
                        value['blob'] if complete else None,pool.queue.archive_case(key,store),
                        None if complete else dict(state=state['state'],error=value.get('error')))
                    if not complete and (stage=='audit' or (state['state']=='error'
                            and 'min() arg is an empty sequence' not in str(value.get('error')))):
                        # A disagreement or invalid audit input stops admission immediately.
                        stopped_on_error=key;exhausted=True
                    del active[key]
                if time.monotonic()>=next_save or not active:
                    ledger.checkpoint();pool.checkpoint()
                    progress=dict(counts=ledger.counts(),active=len(active),drained=exhausted and not active,
                                  stopped_on_error=stopped_on_error)
                    _publish(args.output/'progress.json',progress);print(json.dumps(progress),flush=True)
                    next_save=time.monotonic()+30
                if active:pool.check_workers();time.sleep(.1)
        rows=[]
        for name,stage,state,ref,audit,error in ledger.db.execute('SELECT * FROM results'):
            rows.append(dict(name=name,stage=stage,state=state,
                result=store.load(Blob(**json.loads(ref))) if ref else None,
                audit=json.loads(audit),error=json.loads(error) if error else None))
        _publish(args.output/'results.json',dict(rows=rows,expected_tasks=len(tasks)))
        ledger.db.close()
    with sqlite3.connect(f'file:{source}/ledger.sqlite?mode=ro', uri=True) as db:
        final_rows = db.execute('SELECT * FROM results ORDER BY case_name,stage').fetchall()
    if final_rows != source_rows:
        raise ValueError('Original benchmark stage records changed during audit')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=16)
    run(parser.parse_args())
