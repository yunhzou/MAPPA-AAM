"""Read-only dense bond-count checks of published conversion pilot outputs."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import pickle
import sqlite3
import tempfile
import time

from mappa.mappa_core.core.scheduler import read_blob
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.utils.bonds import bond_event_threshold
from .inputs import load_plan
from .saved_analysis import read_member
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint


def verify_conversion(problem, value, resource, key):
    """Independently count saved branches; never run matching or decoding."""
    import numpy as np
    result = value['decoded']
    assert result['schema'] == 'mappa.decoded_events/v2'
    assert result['event_policy']['schema'] == 'mappa.signed_events/v3'
    assert result['method'] == 'anchored_fragment_conversion'
    catalogue = FinalBranchCatalogue.from_record(problem, result['converted_catalogue'])
    config = result['config']
    edges=[]
    for endpoint in (problem.reactant,problem.product):
        a,b=np.triu_indices(endpoint.atom_count,1)
        thresholds=np.array([bond_event_threshold(endpoint.elements,int(x),int(y),
            default_threshold=config['threshold'],metal_threshold=config['metal_threshold']) for x,y in zip(a,b)])
        edges.append((a,b,endpoint.wbo[a,b],thresholds))
    padded=[np.pad(endpoint.wbo,((0,1),(0,1))) for endpoint in (problem.product,problem.reactant)]
    counts = {}
    reports = {row['family_id']:row for row in result['family_reports']}
    assert len(result['family_reports']) == len(result['converted_branches'])
    assert set(reports) == {b['family_id'] for b in result['converted_branches']}
    for branch in result['converted_branches']:
        family = catalogue.families[branch['family_id']]
        mapping = dict(family.mapping)
        assert len(mapping) == len(set(mapping.values()))
        for pair_id, _ in branch['fragment_paths']:
            pair = result['rematches'][pair_id]
            assert all(mapping[a] == b for a,b in pair['anchors'])
            assert {mapping[a] for a in pair['source_atoms']} <= set(pair['target_atoms'])
        if family.mapping not in counts:
            count = 0
            for size, (a,b,weights,thresholds), other, assignment in zip(
                    (problem.source_atom_count,problem.target_atom_count),edges,padded,
                    (mapping,{p:r for r,p in mapping.items()})):
                # Use the full endpoint dimension, including unmatched atoms.
                images=np.full(size,len(other)-1,dtype=int)
                for atom,image in assignment.items():
                    images[atom]=image
                touched=(images[a]!=len(other)-1)|(images[b]!=len(other)-1)
                count+=np.count_nonzero(touched & (weights-other[images[a],images[b]]>thresholds))
            counts[family.mapping] = int(count)
        report = reports[branch['family_id']]
        assert counts[family.mapping] == report['total'] and report['solver_queries'] == 0
    minimum = min(counts.values(), default=None)
    assert minimum == result['observed_minimum']
    for candidate in result['candidates']:
        assert candidate['total'] <= minimum + config['event_slack']
        assert all(reports[i]['total'] == candidate['total'] for i in candidate['family_ids'])
    return dict(key=key, complete=result['complete'], minimum=minimum,
        branches=len(reports), distinct_mappings=len(counts), candidates=len(result['candidates']),
        internal_events=sum(r['internal_events'] for r in reports.values()),
        incomplete_rematches=sum(not r['complete'] for r in result['rematches']),
        event_histogram=dict(Counter(r['total'] for r in reports.values())),
        timing=value['timing'], resource={k:v for k,v in resource.items() if k!='blob'})

def check(queue_path, output):
    started = time.perf_counter()
    db = sqlite3.connect(f'file:{queue_path.resolve()}?mode=ro&immutable=1', uri=True)
    states = dict(db.execute('SELECT state,count(*) FROM tasks GROUP BY state'))
    if set(states) != {'complete'}:
        raise ValueError(f'Pilot executions are unfinished: {states}')
    summaries, search_outputs = [], 0
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        folder = Path(scratch)
        for key, specification, resource in db.execute('SELECT key,specification,result FROM tasks ORDER BY key'):
            task, resource = _loads_checkpoint(specification), _loads_checkpoint(resource)
            value = read_blob(resource['blob'], scratch)
            if task.operation in ('mappa.mappa_core.core.scheduler.search:search_batch',
                                  'mappa.mappa_core.core.scheduler.search:finalize_cut'):
                assert value['graph'].contexts
                search_outputs += 1
                continue
            if task.operation not in ('rebenchmark.decoding:decode_seed','rebenchmark.operations:convert_seed'):
                raise ValueError(f'Unexpected pilot operation: {task.operation}')
            assert key == f'{value["task_id"]:05d}/cut/{value["cut"]:05d}/seed/{value["seed"]-1:02d}'
            if task.operation == 'rebenchmark.operations:convert_seed':
                from mappa.mappa_core.utils.io.records import problem_from_record
                problem = problem_from_record(read_blob(value['input'],scratch))
                assert value['search_rerun']
                assert value['finalized_sha256'] == value['finalized'].sha256
                read_blob(value['finalized'],scratch)
            else:
                for name, member in value['source'].items():
                    (folder/name).write_bytes(read_member(member))
                problem = load_plan(folder).problem
                assert not value['search_rerun']
            summaries.append(verify_conversion(problem, value, resource, key))
    db.close()
    summary = dict(search_outputs_verified=search_outputs, tasks=len(summaries), complete=sum(r['complete'] for r in summaries),
        branches=sum(r['branches'] for r in summaries),
        internal_events=sum(r['internal_events'] for r in summaries),
        seconds=time.perf_counter()-started, rows=summaries, verified=True,
        scope='All saved pilot branches independently counted; anchors, atom sets and publication hashes verified. No search or decoding repeated.')
    output.write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    check(args.queue, args.output)
