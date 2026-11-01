"""Prepare disjoint input batches for a fresh run or a stopped campaign."""
import argparse
from collections import defaultdict
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sqlite3

from mappa import AAMSearchConfig


def prepare_fresh(index_paths, output, workers):
    """Balance immutable inputs; import no results, attempts or cache snapshots."""
    if not workers or any(count <= 0 for count in workers):
        raise ValueError('Each batch requires a positive worker count')
    sources = [Path(path) for path in index_paths]
    indexes = [json.loads(path.read_bytes()) for path in sources]
    settings = [{key: value for key, value in index.items() if key != 'directions'}
                for index in indexes]
    if not settings or any(item != settings[0] for item in settings):
        raise ValueError('Input manifests must use identical scientific settings')
    metadata = settings[0]
    config = AAMSearchConfig(**metadata['config'])
    metadata = dict(metadata, config=asdict(config))
    reactions = defaultdict(list)
    identities = set()
    for index in indexes:
        for row in index['directions']:
            identity = row['case'], row['reversed']
            if identity in identities:
                raise ValueError(f'Duplicate reaction/direction input: {identity}')
            identities.add(identity)
            reactions[row['case']].append(row)
    if any(len(rows) != 2 or {row['reversed'] for row in rows} != {False, True}
           for rows in reactions.values()):
        raise ValueError('Every reaction must have both directions')
    def size(case):
        return sum((row['input'] if 'input' in row else row['blob'])['size']
                   for row in reactions[case])
    batches = [[] for _ in workers]
    loads = [0] * len(workers)
    for case in sorted(reactions, key=lambda case: (-size(case), case)):
        batch = min(range(len(workers)), key=lambda i: (loads[i]/workers[i], i))
        batches[batch].extend(sorted(reactions[case], key=lambda row: row['reversed']))
        loads[batch] += size(case)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    summaries = []
    for number, rows in enumerate(batches):
        folder = output/str(number)
        folder.mkdir()
        raw = json.dumps(dict(metadata, directions=rows)).encode()
        (folder/'index.json').write_bytes(raw)
        summaries.append(dict(batch=number, workers=workers[number], cases=len(rows)//2,
            directions=len(rows), input_bytes=loads[number],
            index_sha256=hashlib.sha256(raw).hexdigest()))
    manifest = dict(cases=len(reactions), directions=len(identities), config=asdict(config),
        batches=summaries, fresh=True, reused_scientific_results=0,
        ordering='descending prepared-input bytes; greedy load per planned worker',
        sources=[dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                 for path in sources])
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


def prepare(index_path, previous, output, count=40, additional_ledgers=()):
    index = json.loads(Path(index_path).read_text())
    previous = Path(previous)
    plan = json.loads((previous/'plan.json').read_text())
    assert plan['input_index_sha256'] == hashlib.sha256(Path(index_path).read_bytes()).hexdigest()
    with sqlite3.connect(f'file:{previous}/ledger.sqlite?mode=ro', uri=True) as source:
        records = source.execute('SELECT * FROM results').fetchall()
        schema = source.execute("SELECT sql FROM sqlite_master WHERE name='results'").fetchone()[0]
    merged = {(row[0], row[1]): row for row in records}
    for ledger in additional_ledgers:
        with sqlite3.connect(f'file:{ledger}?mode=ro', uri=True) as source:
            for row in source.execute('SELECT * FROM results'):
                key = row[0], row[1]
                if key in merged and merged[key] != row:
                    raise ValueError(f'Conflicting checkpoint results for {key}')
                merged[key] = row
    records = list(merged.values())
    for batch in range(count):
        root = Path(output)/str(batch)
        root.mkdir(parents=True, exist_ok=False)
        rows = [row for row in index['directions'] if row['case'] % count == batch]
        # A batch manifest contains only its own scientific input rows.
        subset = dict(config=index['config'], directions=rows)
        path = root/'index.json'
        path.write_text(json.dumps(subset))
        result = root/'result'; result.mkdir()
        names = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}" for row in rows}
        updated = dict(plan, input_index=str(path.resolve()),
                       input_index_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       cases=len({row['case'] for row in rows}), directions=len(rows))
        (result/'plan.json').write_text(json.dumps(updated,indent=2)+'\n')
        with sqlite3.connect(result/'ledger.sqlite') as destination:
            destination.execute(schema)
            destination.executemany('INSERT INTO results VALUES (?,?,?,?,?,?)',
                                    [record for record in records if record[0] in names])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='operation', required=True)
    fresh = commands.add_parser('fresh', help='New manifests without any saved scientific results')
    fresh.add_argument('output'); fresh.add_argument('indexes', nargs='+')
    fresh.add_argument('--workers', type=int, nargs='+', required=True)
    resume = commands.add_parser('resume', help='Repartition completed results from a stopped campaign')
    resume.add_argument('index'); resume.add_argument('previous'); resume.add_argument('output')
    resume.add_argument('--count',type=int,default=40)
    resume.add_argument('--additional-ledger',action='append',default=[])
    args = parser.parse_args()
    if args.operation == 'fresh':
        result = prepare_fresh(args.indexes, args.output, args.workers)
        print(json.dumps({key: value for key, value in result.items() if key != 'sources'}))
    else:
        prepare(args.index, args.previous, args.output,count=args.count,
                additional_ledgers=args.additional_ledger)
