"""Prepare disjoint input batches for a fresh run or a stopped campaign."""
import argparse
from collections import defaultdict
from contextlib import ExitStack, closing
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sqlite3

from mappa import AAMSearchConfig
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blobs
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from enzymemap import decoding_policy


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
    metadata = dict(metadata, config=asdict(config), decode_config=asdict(decoding_policy(metadata)))
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
        decode_config=metadata['decode_config'],
        batches=summaries, fresh=True, reused_scientific_results=0,
        ordering='descending prepared-input bytes; greedy load per planned worker',
        sources=[dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                 for path in sources])
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


def _direction_rows(index):
    rows = {}
    for row in index['directions']:
        if not isinstance(row['case'], int) or type(row['reversed']) is not bool:
            raise ValueError('Directions require integer cases and boolean reversed flags')
        name = f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
        if name in rows:
            raise ValueError(f'Duplicate reaction/direction input: {name}')
        rows[name] = row
    return rows


def _decoding_catalogues(source, owners):
    """Read a complete frozen campaign and its small archived attempt records."""
    from enzymemap import recovery_owner

    raw = (source/'manifest.json').read_bytes()
    manifest = json.loads(raw)
    launch_raw = (source/'launch.json').read_bytes()
    launch = json.loads(launch_raw)
    config = asdict(AAMSearchConfig(**manifest['config']))
    policy = decoding_policy(manifest)
    expected, source_indexes = {}, []
    for origin in manifest['sources']:
        path = Path(origin['path'])
        index_raw = path.read_bytes()
        digest = hashlib.sha256(index_raw).hexdigest()
        if digest != origin['sha256']:
            raise ValueError('Original source index checksum differs')
        index = json.loads(index_raw)
        if asdict(AAMSearchConfig(**index['config'])) != config or decoding_policy(index) != policy:
            raise ValueError('Original source configuration differs')
        rows = _direction_rows(index)
        if expected.keys() & rows.keys():
            raise ValueError('Duplicate reaction/direction across original source indexes')
        expected.update(rows)
        source_indexes.append(dict(path=str(path.resolve()), sha256=digest))
    cases = {row['case'] for row in expected.values()}
    if (len(expected) != manifest['directions'] or len(cases) != manifest['cases']
            or len(expected) != 2*len(cases)):
        raise ValueError('Original campaign must contain both directions of every reaction')
    selected, memory_failures, origins, identities, batch_ids = [], set(), [], set(), set()
    for batch in manifest['batches']:
        number = batch['batch']
        if not isinstance(number, int) or number < 0 or number in batch_ids:
            raise ValueError('Original campaign batch identities must be unique nonnegative integers')
        batch_ids.add(number)
        directory = source/str(number)
        result = directory/'result'
        recovery_owner(result, owners)
        index_raw = (directory/'index.json').read_bytes()
        plan_raw = (result/'plan.json').read_bytes()
        progress_raw = (result/'progress.json').read_bytes()
        index, plan, progress = map(json.loads, (index_raw, plan_raw, progress_raw))
        digest = hashlib.sha256(index_raw).hexdigest()
        if (digest != batch['index_sha256'] or digest != plan['input_index_sha256']
                or index['config'] != plan['config']
                or asdict(AAMSearchConfig(**index['config'])) != config
                or decoding_policy(index) != policy or decoding_policy(plan) != policy):
            raise ValueError('Original campaign input or configuration differs')
        if progress['source'] != launch['source_identity']:
            raise ValueError('Original campaign source identities differ')
        rows = _direction_rows(index)
        if identities.intersection(rows):
            raise ValueError('Duplicate reaction/direction across original batches')
        identities.update(rows)
        if any(expected.get(name) != row for name, row in rows.items()):
            raise ValueError('Original batch input differs from its source index')
        count = len({row['case'] for row in rows.values()})
        if any(value['directions'] != len(rows) or value['cases'] != count
               for value in (batch, plan)):
            raise ValueError('Original batch direction or reaction counts differ')
        path = (result/'ledger.sqlite').resolve()
        with closing(sqlite3.connect(f'{path.as_uri()}?mode=ro', uri=True)) as ledger:
            records = {(row[0], row[1]): row for row in ledger.execute('SELECT * FROM results')}
            retries = ledger.execute('SELECT case_name,stage,state,reference,audit,error FROM retries').fetchall()
        if any(name not in rows for name, _ in records):
            raise ValueError('Original ledger contains directions outside its batch')
        for name, row in rows.items():
            for stage in ('matching', 'deduplication'):
                record = records.get((name, stage))
                if record is None or record[2] != 'complete' or not record[3]:
                    raise ValueError('Fresh decoding requires complete matching and reduction')
            catalogue = Blob(**json.loads(records[name, 'deduplication'][3]))
            selected.append((name, dict(row, catalogue=asdict(catalogue))))
        # A successful retry can still have earlier memory-limit attempts.
        # Read archived audits only; catalogue payloads remain unopened.
        audits = {}
        for name, stage, state, reference, audit, error in [*records.values(), *retries]:
            if name not in rows:
                raise ValueError('Original retry belongs to another batch')
            if stage != 'decoding':
                continue
            if state not in ('complete', 'error') or not audit:
                raise ValueError('Original decoding attempts must be terminal and archived')
            blob = Blob(**json.loads(audit))
            audits.setdefault(blob, set()).add(name)
        for blob, audit in read_blobs(sorted(audits, key=lambda ref: (ref.path, ref.offset)), None):
            if (any(task[4] == 'memory_limit' for task in audit['tasks'])
                    or any(attempt[5] == 'memory_limit' for attempt in audit['attempts'])):
                memory_failures.update(audits[blob])
        origins.append(dict(batch=number, input_index_sha256=digest,
            plan_sha256=hashlib.sha256(plan_raw).hexdigest(),
            ledger_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            queue_sha256=hashlib.sha256((result/'bundles/queue.sqlite').read_bytes()).hexdigest(),
            progress_sha256=hashlib.sha256(progress_raw).hexdigest()))
    if identities != expected.keys():
        raise ValueError('Original batches do not cover the complete source direction set')
    provenance = dict(campaign=str(source), manifest_sha256=hashlib.sha256(raw).hexdigest(),
        launch_sha256=hashlib.sha256(launch_raw).hexdigest(), source_commit=launch['source_commit'],
        source_identity=launch['source_identity'], source_indexes=source_indexes,
        original_batches=origins)
    return config, policy, selected, memory_failures, provenance


def prepare_decoding(source, output, workers, *, memory_failures_workers=(), decode_config=None):
    """Prepare a fresh decode pass from every saved catalogue, importing no result rows."""
    if not workers or any(count <= 0 for count in (*workers, *memory_failures_workers)):
        raise ValueError('Each batch requires a positive worker count')
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError('Fresh decoding requires a separate output directory')
    if output.exists():
        raise FileExistsError(output)
    with ExitStack() as owners:
        config, original_policy, selected, failures, provenance = _decoding_catalogues(source, owners)
        policy = decode_config or original_policy
        if not isinstance(policy, EventDecodeConfig):
            raise TypeError('decode_config must be EventDecodeConfig')
        output.mkdir(parents=True)
        metadata = dict(config=config, decode_config=asdict(policy), decoding_source=provenance)
        groups = [('ordinary', workers, selected)]
        if memory_failures_workers:
            groups = [('ordinary', workers, [(name, row) for name, row in selected if name not in failures]),
                ('memory_failures', memory_failures_workers,
                 [(name, row) for name, row in selected if name in failures])]
        batches, summaries = [], []
        for group, counts, rows in groups:
            buckets, loads = [[] for _ in counts], [0]*len(counts)
            for name, row in sorted(rows, key=lambda item: (-item[1]['catalogue']['size'], item[0])):
                number = min(range(len(counts)), key=lambda i: (loads[i]/counts[i], i))
                buckets[number].append(row)
                loads[number] += row['catalogue']['size']
            for count, bucket, size in zip(counts, buckets, loads):
                raw = json.dumps(dict(metadata, directions=bucket)).encode()
                summaries.append(dict(batch=len(batches), workers=count, resource_group=group,
                    cases=len({row['case'] for row in bucket}), directions=len(bucket),
                    catalogue_bytes=size, index_sha256=hashlib.sha256(raw).hexdigest()))
                batches.append(raw)
        for number, raw in enumerate(batches):
            folder = output/str(number)
            folder.mkdir()
            (folder/'index.json').write_bytes(raw)
        manifest = dict(metadata, cases=len({row['case'] for _, row in selected}),
            directions=len(selected), batches=summaries, fresh=True, operation='decoding',
            reused_catalogues=len(selected), imported_stage_results=0,
            prior_memory_failure_directions=sorted(failures),
            ordering='descending catalogue bytes; greedy load per planned worker within each resource group')
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
        subset = dict(index, directions=rows)
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
    decoding = commands.add_parser('decoding', help='Fresh decode manifests from a completed matching campaign')
    decoding.add_argument('source'); decoding.add_argument('output')
    decoding.add_argument('--workers', type=int, nargs='+', required=True)
    decoding.add_argument('--memory-failures-workers', type=int, nargs='+', default=[])
    resume = commands.add_parser('resume', help='Repartition completed results from a stopped campaign')
    resume.add_argument('index'); resume.add_argument('previous'); resume.add_argument('output')
    resume.add_argument('--count',type=int,default=40)
    resume.add_argument('--additional-ledger',action='append',default=[])
    args = parser.parse_args()
    if args.operation == 'fresh':
        result = prepare_fresh(args.indexes, args.output, args.workers)
        print(json.dumps({key: value for key, value in result.items() if key != 'sources'}))
    elif args.operation == 'decoding':
        result = prepare_decoding(args.source, args.output, args.workers,
                                  memory_failures_workers=args.memory_failures_workers)
        print(json.dumps({key: value for key, value in result.items()
                          if key not in ('decoding_source', 'prior_memory_failure_directions')}))
    else:
        prepare(args.index, args.previous, args.output,count=args.count,
                additional_ledgers=args.additional_ledger)
