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
from mappa.final_branches import FinalBranchCatalogue
from mappa.frag import bond_event_threshold
from .campaign import load_plan
from .saved_analysis import read_member
from mappa.artifacts import _loads_checkpoint


def check(queue_path, output):
    started = time.perf_counter()
    db = sqlite3.connect(f'file:{queue_path.resolve()}?mode=ro&immutable=1', uri=True)
    states = dict(db.execute('SELECT state,count(*) FROM tasks GROUP BY state'))
    if set(states) != {'complete'}:
        raise ValueError(f'Pilot executions are unfinished: {states}')
    summaries = []
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        folder = Path(scratch)
        for key, specification, resource in db.execute('SELECT key,specification,result FROM tasks ORDER BY key'):
            task, resource = _loads_checkpoint(specification), _loads_checkpoint(resource)
            value = read_blob(resource['blob'], scratch)
            assert task.operation == 'rebenchmark.decoding:decode_seed'
            assert key == f'{value["task_id"]:05d}/cut/{value["cut"]:05d}/seed/{value["seed"]-1:02d}'
            for name, member in value['source'].items():
                (folder/name).write_bytes(read_member(member))
            problem = load_plan(folder).problem
            result = value['decoded']
            assert result['schema'] == 'mappa.decoded_events/v2'
            assert result['method'] == 'anchored_fragment_conversion' and not value['search_rerun']
            catalogue = FinalBranchCatalogue.from_record(problem, result['converted_catalogue'])
            config = result['config']
            pairs = [tuple((a, b, ep.wbo[a,b], bond_event_threshold(ep.elements, a, b,
                default_threshold=config['threshold'], metal_threshold=config['metal_threshold']))
                for a in range(ep.atom_count) for b in range(a+1, ep.atom_count))
                for ep in (problem.reactant, problem.product)]
            counts = {}
            reports = {row['family_id']:row for row in result['family_reports']}
            assert len(reports) == len(catalogue.families) == len(result['converted_branches'])
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
                    for edges, other, assignment in ((pairs[0], problem.product, mapping),
                            (pairs[1], problem.reactant, {p:r for r,p in mapping.items()})):
                        for a,b,weight,threshold in edges:
                            if a not in assignment and b not in assignment:
                                continue
                            opposite = other.wbo[assignment[a],assignment[b]] if a in assignment and b in assignment else 0.
                            count += weight - opposite >= threshold
                    counts[family.mapping] = int(count)
                report = reports[branch['family_id']]
                assert counts[family.mapping] == report['total'] and report['solver_queries'] == 0
            minimum = min(counts.values(), default=None)
            assert minimum == result['observed_minimum']
            for candidate in result['candidates']:
                assert candidate['total'] <= minimum + config['event_slack']
                assert all(reports[i]['total'] == candidate['total'] for i in candidate['family_ids'])
            summaries.append(dict(key=key, complete=result['complete'], minimum=minimum,
                branches=len(reports), distinct_mappings=len(counts), candidates=len(result['candidates']),
                internal_events=sum(r['internal_events'] for r in reports.values()),
                incomplete_rematches=sum(not r['complete'] for r in result['rematches']),
                event_histogram=dict(Counter(r['total'] for r in reports.values())),
                timing=value['timing'], resource={k:v for k,v in resource.items() if k!='blob'}))
    db.close()
    summary = dict(tasks=len(summaries), complete=sum(r['complete'] for r in summaries),
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
