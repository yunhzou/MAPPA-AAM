"""Replay saved fragment-response fixtures without rerunning atom matching."""
import argparse
from dataclasses import asdict
import gzip
import json
from pathlib import Path
import pickle
from statistics import median
import sys
from time import perf_counter, process_time

from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.core.symmetry.minimum import FragmentSymmetryEvents
from mappa.mappa_core.core.symmetry.objective import FragmentObjective


def timed(operation):
    wall, cpu = perf_counter(), process_time()
    result = operation()
    return result, dict(wall_seconds=perf_counter()-wall, cpu_seconds=process_time()-cpu)


def measure(fragment, *, repeats=5):
    timings = []
    for _ in range(repeats):
        response, elapsed = timed(lambda: FragmentSymmetryEvents(fragment))
        timings.append(elapsed)
    response.deadline = perf_counter() + 5.
    minimum = response.solve()
    response.deadline = None
    views = (response.select(minimum.event_count), response.select(minimum.event_count+1),
             response.select(minimum.event_count, minimum.event_count+2))
    payload = pickle.dumps(views, protocol=5)
    restored = pickle.loads(payload)
    assert all(view.response is restored[0].response for view in restored)
    queries = []
    objective = FragmentObjective(fragment)
    for view in restored:
        samples = []
        for _ in range(repeats):
            answer, elapsed = timed(view.find)
            samples.append(elapsed)
        if answer.status == 'found':
            actual = objective.count(objective.images(answer.witness))
            assert actual == answer.event_count
        queries.append(dict(bounds=(view.lower,view.upper), status=answer.status,
            event_count=answer.event_count, witness=answer.witness,
            median={key:median(row[key] for row in samples) for key in elapsed}))
    return dict(atoms=len(fragment.source_atoms), minimum=asdict(minimum),
        preparation_median={key:median(row[key] for row in timings) for key in elapsed},
        fixed_events=response.fixed_events, variable_terms=len(response.pairs),
        structural_operations=response.group_order, queries=queries,
        compressed_bytes=len(gzip.compress(payload,mtime=0)))


def compare(fragment, reference, *, repeats=5):
    """Interleave the accepted standalone beta and production preparation."""
    positions = {atom:i for i,atom in enumerate(fragment.target_atoms)}
    source_positions = {atom:i for i,atom in enumerate(fragment.source_atoms)}
    source, target = fragment.source_symmetry, fragment.target_symmetry
    old = reference.Fragment(source.elements,target.elements,source.weights,target.weights,
        tuple(positions[p] for _,p in fragment.witness),
        tuple((source_positions[a],positions[b]) for a,b in fragment.relation.anchors))
    operations = dict(beta=lambda:reference.PreparedFragment(old),
                      mappa=lambda:FragmentSymmetryEvents(fragment))
    samples = {name:[] for name in operations}
    for operation in operations.values():
        operation()
    for repeat in range(repeats):
        for name in (tuple(operations) if repeat % 2 == 0 else tuple(reversed(operations))):
            response, elapsed = timed(operations[name])
            samples[name].append(elapsed)
    return {name:{key:median(row[key] for row in rows) for key in rows[0]}
            for name,rows in samples.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--reference',type=Path,help='Accepted standalone symmetry-event-lab checkout')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    data = json.loads(gzip.decompress(args.inputs.read_bytes()))
    reference = None
    if args.reference is not None:
        sys.path.insert(0,str(args.reference.resolve()))
        import prototype as reference
    rows = []
    for saved in data['fragments']:
        fragment = Fragment.from_record(saved['fragment'])
        result = measure(fragment,repeats=args.repeats)
        if reference is not None:
            result['interleaved_preparation'] = compare(fragment,reference,repeats=args.repeats)
        rows.append({**saved,'measurement':result})
        print(saved['name'],result['preparation_median'],flush=True)
    args.output.write_bytes(gzip.compress(json.dumps(dict(
        scope=data['scope'], repeats=args.repeats, fragments=rows),separators=(',',':')).encode(),mtime=0))


if __name__ == '__main__':
    main()
