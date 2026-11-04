"""Fresh EnzymeMap passes over prepared inputs, using one persistent MAPPA pool."""
import argparse
from contextlib import ExitStack, closing
from dataclasses import asdict
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import sqlite3
import tempfile
import time
import uuid

from mappa import AAMSearchConfig, EventDecodeConfig, query_path
from mappa.execution import BundleStore, StageResults, WorkPool
from mappa.execution.identity import code_identity
from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
from mappa.mappa_core.utils.io.records import problem_from_record


def decoding_policy(index):
    """Read the declared event policy, using defaults only for omitted settings."""
    return EventDecodeConfig(**index.get('decode_config', dict(
        threshold=index['config']['event_threshold'],
        metal_threshold=index['config']['metal_event_threshold'])))


def prepared_problem(row):
    if 'input' in row:
        return problem_from_record(read_blob(Blob(**row['input']), None))
    # The preparation index points to a saved task list. Only its endpoint-input
    # reference is read; none of those retired task operations is executed.
    preparation = read_blob(Blob(**row['blob']), None)
    source = read_blob(preparation[0].payload['input'], None)
    return problem_from_record(source)


def verify(payload, workspace):
    saved = workspace.load(payload['catalogue'])
    reference = read_blob(payload['reference'], workspace.scratch)
    target = dict(reference['mapping'])
    if payload['reversed']:
        target = {product: reactant for reactant, product in target.items()}
    unknown = False
    for fid, family in enumerate(saved.catalogue.families):
        status, evidence = query_path(family.as_path(saved.catalogue.problem), saved.catalogue.problem,
            target, source_atoms=tuple(sorted(target)), complete_reference=True, timeout_ms=1000)
        if status == 'recovered':
            return dict(status=status, family=fid, evidence=evidence, heavy_atoms=len(target),
                        policy='exact annotated heavy-atom indices; no endpoint-symmetry relaxation')
        unknown |= status == 'unknown'
    return dict(status='unknown' if unknown else 'not_recovered', heavy_atoms=len(target),
                policy='exact annotated heavy-atom indices; no endpoint-symmetry relaxation')


class Ledger:
    """Local transactions; one published SQLite snapshot for all reactions/stages."""
    def __init__(self, output, scratch, *, resume):
        self.output = Path(output)
        self.path = Path(scratch)/'ledger.sqlite'
        saved = self.output/'ledger.sqlite'
        if resume:
            shutil.copyfile(saved, self.path)
        elif saved.exists():
            raise ValueError('Existing results require --resume')
        self.db = sqlite3.connect(self.path)
        self.db.execute('CREATE TABLE IF NOT EXISTS results (case_name TEXT, stage TEXT, state TEXT, reference TEXT, audit TEXT, error TEXT, PRIMARY KEY(case_name,stage))')
        self.db.execute('CREATE TABLE IF NOT EXISTS retries (utc REAL, reason TEXT, case_name TEXT, stage TEXT, state TEXT, reference TEXT, audit TEXT, error TEXT)')
        self.db.commit()

    def get(self, case, stage):
        return self.db.execute('SELECT state,reference FROM results WHERE case_name=? AND stage=?', (case,stage)).fetchone()

    def put(self, case, stage, state, reference=None, audit=None, error=None):
        with self.db:
            self.db.execute('INSERT INTO results VALUES (?,?,?,?,?,?)',
                (case, stage, state, json.dumps(asdict(reference)) if reference else None,
                 json.dumps(asdict(audit)) if audit else None, json.dumps(error) if error else None))

    def retry_audit(self, case, stage):
        row = self.db.execute('SELECT audit FROM results WHERE case_name=? AND stage=? AND state=?',
                              (case, stage, 'error')).fetchone()
        return Blob(**json.loads(row[0])) if row and row[0] else None

    def retire_failure(self, case, stage, reason):
        with self.db:
            self.db.execute('INSERT INTO retries SELECT ?,?,* FROM results WHERE case_name=? AND stage=? AND state=?',
                            (time.time(), reason, case, stage, 'error'))
            self.db.execute('DELETE FROM results WHERE case_name=? AND stage=? AND state=?', (case, stage, 'error'))

    def checkpoint(self):
        local = self.path.with_name('snapshot.sqlite')
        with sqlite3.connect(local) as destination:
            self.db.backup(destination)
        pending = self.output/'ledger.pending.sqlite'
        shutil.copyfile(local, pending)
        pending.replace(self.output/'ledger.sqlite')

    def counts(self):
        return [dict(stage=stage,state=state,count=count) for stage,state,count in
                self.db.execute('SELECT stage,state,count(*) FROM results GROUP BY stage,state')]


def require_drained(output):
    queue_path = (Path(output)/'bundles/queue.sqlite').resolve()
    with closing(sqlite3.connect(f'{queue_path.as_uri()}?mode=ro', uri=True)) as queue:
        if queue.execute('SELECT EXISTS(SELECT 1 FROM tasks) OR EXISTS(SELECT 1 FROM attempts)').fetchone()[0]:
            raise ValueError('Execution requires a drained saved queue; unarchived tasks or attempts remain')


def recovery_owner(output, owners):
    owner = owners.enter_context((Path(output)/'bundles/cache/.owner.lock').open('rb'))
    try:
        fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError('Original campaign or predecessor has an active coordinator or decoding recovery') from None
    require_drained(output)


def recovery_registration(output, plan):
    return dict(schema='mappa.decoding-recovery-owner/v1', output=str(Path(output).resolve()),
                plan_sha256=hashlib.sha256(json.dumps(plan, sort_keys=True).encode()).hexdigest())


def campaign_decoding_recovery_plan(source, output, memory_bytes, owners):
    """Select terminal decode failures without loading catalogues or changing their campaign."""
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError('Decoding recovery requires a separate output directory')
    raw = (source/'manifest.json').read_bytes()
    manifest = json.loads(raw)
    decoding_config = asdict(decoding_policy(manifest))
    selected, origins, identities, task_keys = [], [], set(), set()
    for batch in manifest['batches']:
        directory = source/str(batch['batch'])
        result = directory/'result'
        # Hold original ownership through recovery: another output must not
        # retry these same tasks or resume their source campaign concurrently.
        recovery_owner(result, owners)
        index_raw = (directory/'index.json').read_bytes()
        index = json.loads(index_raw)
        plan_raw = (result/'plan.json').read_bytes()
        plan = json.loads(plan_raw)
        digest = hashlib.sha256(index_raw).hexdigest()
        if (digest != batch['index_sha256'] or digest != plan['input_index_sha256']
                or index['config'] != plan['config'] or index['config'] != manifest['config']
                or asdict(decoding_policy(index)) != decoding_config
                or asdict(decoding_policy(plan)) != decoding_config):
            raise ValueError('Original campaign input or configuration differs')
        names = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}": row
                 for row in index['directions']}
        if len(names) != len(index['directions']) or identities.intersection(names):
            raise ValueError('Original campaign contains duplicate reaction/direction inputs')
        identities.update(names)
        path = (result/'ledger.sqlite').resolve()
        with closing(sqlite3.connect(f'{path.as_uri()}?mode=ro', uri=True)) as ledger:
            records = {(row[0], row[1]): row for row in ledger.execute('SELECT * FROM results')}
        stages = ('matching', 'deduplication', 'decoding', 'verification')
        if (set(records) != {(name, stage) for name in names for stage in stages}
                or any(row[2] not in ('complete', 'error') for row in records.values())):
            raise ValueError('Original campaign has incomplete stage coverage')
        for name, row in names.items():
            if any(records[name, stage][2] != 'complete' for stage in ('matching', 'deduplication')):
                raise ValueError('Decoding recovery requires completed matching and reduction')
            saved = records[name, 'decoding']
            if saved[2] == 'complete':
                continue
            audit = Blob(**json.loads(saved[4]))
            record = read_blob(audit, None)
            if len(record['tasks']) != 1 or record['tasks'][0][4] != 'memory_limit':
                raise ValueError('Decoding recovery selects only archived memory failures')
            task_row = record['tasks'][0]
            task = _loads_checkpoint(task_row[1])
            catalogue = Blob(**json.loads(records[name, 'deduplication'][3]))
            if (task.operation != 'mappa.mappa_core.core.scheduler.decoding:decode_catalogue'
                    or task.payload['catalogue'] != catalogue):
                raise ValueError('Original decoding task differs from its saved catalogue')
            if asdict(task.payload['config']) != decoding_config:
                raise ValueError('Original decoding configuration differs from the campaign configuration')
            if task_row[3] != memory_bytes or task.memory_bytes != memory_bytes:
                raise ValueError('Decoding recovery must preserve the original memory guard')
            if record['case'][3] >= record['case'][1]:
                raise ValueError('Decoding recovery has no remaining case budget')
            if task.key in task_keys:
                raise ValueError('Original decoding task appears more than once')
            task_keys.add(task.key)
            selected.append(dict(case=row['case'], reversed=row['reversed'], batch=batch['batch'],
                audit=asdict(audit), catalogue=asdict(catalogue), task=task.key,
                budget_seconds=record['case'][1], spent_seconds=record['case'][3],
                decoding_config=asdict(task.payload['config'])))
        origins.append(dict(batch=batch['batch'], input_index_sha256=digest,
            plan_sha256=hashlib.sha256(plan_raw).hexdigest(),
            ledger_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            source=json.loads((result/'progress.json').read_text())['source']))
    if len(identities) != manifest['directions'] or not selected:
        raise ValueError('Original campaign scope differs or contains no decoding memory failures')
    original_sources = {origin['source'] for origin in origins}
    if 'code_identity' in manifest:
        original_sources.add(manifest['code_identity'])
    launch = source/'launch.json'
    if launch.exists():
        original_sources.add(json.loads(launch.read_text())['code_identity'])
    if len(original_sources) != 1:
        raise ValueError('Original campaign source identities differ')
    requested = dict(schema='mappa.decoding-recovery/v1', original_campaign=str(source),
        manifest_sha256=hashlib.sha256(raw).hexdigest(), original_batches=origins,
        original_source=original_sources.pop(),
        decoder_source=code_identity(), memory_bytes=memory_bytes, directions=len(selected),
        cases=len({row['case'] for row in selected}), directions_to_recover=selected)
    return selected, requested


def decoding_recovery_plan(source, output, memory_bytes, owners):
    """Validate a registered recovery chain and continue only its latest failures."""
    source, output = Path(source).resolve(), Path(output).resolve()
    predecessors = []
    directory = source
    while not (directory/'manifest.json').exists():
        plan = json.loads((directory/'plan.json').read_text())
        if plan.get('schema') != 'mappa.decoding-recovery/v1':
            raise ValueError('Recovery source must be an original campaign or a decoding recovery')
        if directory in predecessors:
            raise ValueError('Decoding recovery lineage contains a cycle')
        predecessors.append(directory)
        if 'recovery_source' in plan:
            if plan.get('source_kind') != 'decoding_recovery':
                raise ValueError('Decoding recovery source kind differs')
            directory = Path(plan['recovery_source']).resolve()
        else:
            directory = Path(plan['original_campaign']).resolve()
    selected, requested = campaign_decoding_recovery_plan(directory, output, memory_bytes, owners)
    registration = directory/'decoding-recovery.json'
    for directory in reversed(predecessors):
        if output.is_relative_to(directory) or directory.is_relative_to(output):
            raise ValueError('Decoding recovery requires a separate output directory')
        recovery_owner(directory, owners)
        plan_raw = (directory/'plan.json').read_bytes()
        plan = json.loads(plan_raw)
        if json.loads(registration.read_text()) != recovery_registration(directory, plan):
            raise ValueError('Predecessor differs from its registered decoding recovery')
        expected = dict(requested, decoder_source=plan['decoder_source'])
        if plan != expected:
            raise ValueError('Predecessor input, configuration or recovery lineage differs')
        if json.loads((directory/'progress.json').read_text())['source'] != plan['decoder_source']:
            raise ValueError('Predecessor decoder source identity differs')
        ledger_path = (directory/'ledger.sqlite').resolve()
        with closing(sqlite3.connect(f'{ledger_path.as_uri()}?mode=ro', uri=True)) as ledger:
            records = {(row[0], row[1]): row for row in ledger.execute('SELECT * FROM results')}
        names = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}": row
                 for row in selected}
        if (set(records) != {(name, 'decoding') for name in names}
                or any(row[2] not in ('complete', 'error') for row in records.values())):
            raise ValueError('Predecessor recovery has incomplete terminal scope')
        remaining = []
        for name, row in names.items():
            saved = records[name, 'decoding']
            previous = read_blob(Blob(**row['audit']), None)
            audit = Blob(**json.loads(saved[4]))
            current = read_blob(audit, None)
            old_task = previous['tasks'][0]
            old_attempts = {attempt[0]: attempt for attempt in previous['attempts']}
            attempts = {attempt[0]: attempt for attempt in current['attempts']}
            if (len(current['tasks']) != 1 or current['tasks'][0][:4] != old_task[:4]
                    or current['case'][:2] != previous['case'][:2]
                    or current['case'][3] < previous['case'][3]
                    or current['case'][2] is not None or current['case'][4] is not None
                    or current['dependencies'] != previous['dependencies']
                    or current['retries'] != previous['retries']
                    or len(old_attempts) != len(previous['attempts'])
                    or len(attempts) != len(current['attempts'])
                    or len(attempts) != len(old_attempts)+1
                    or any(attempts.get(token) != attempt for token, attempt in old_attempts.items())):
                raise ValueError('Predecessor changed the original task, attempts or consumed budget')
            task = current['tasks'][0]
            token, = attempts.keys()-old_attempts.keys()
            attempt = attempts[token]
            if (attempt[1] != task[0] or attempt[2] != task[5] or token != task[6]
                    or attempt[4] is None or attempt[5] != task[4] or attempt[6] != task[8]):
                raise ValueError('Predecessor task differs from its terminal attempt')
            state = task[4]
            if saved[2] == 'complete':
                result = _loads_checkpoint(task[8])
                if (state != 'complete' or saved[3] is None
                        or result['blob'] != Blob(**json.loads(saved[3]))):
                    raise ValueError('Predecessor successful decoding lacks its completed result')
                continue
            if state != 'memory_limit':
                raise ValueError('Decoding recovery selects only archived memory failures')
            if current['case'][3] >= current['case'][1]:
                raise ValueError('Decoding recovery has no remaining case budget')
            remaining.append(dict(row, audit=asdict(audit), spent_seconds=current['case'][3]))
        if not remaining:
            raise ValueError('Predecessor contains no decoding memory failures')
        history = [*plan.get('predecessors', []), dict(path=str(directory),
            plan_sha256=hashlib.sha256(plan_raw).hexdigest(),
            ledger_sha256=hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
            decoder_source=plan['decoder_source'])]
        selected = remaining
        requested = dict(plan, source_kind='decoding_recovery', recovery_source=str(directory),
            predecessors=history, decoder_source=code_identity(), directions=len(selected),
            cases=len({row['case'] for row in selected}), directions_to_recover=selected)
        registration = directory/'decoding-recovery.json'
    return selected, requested


def run(index_path, output, *, workers=32, window=128, memory_bytes=4*2**30,
        case_seconds=3600., local=False, limit=None, resume=False, run_seconds=12600., stop_after=None,
        retry_memory_limits=False, verify_only=False, recover_decoding=False, decode_only=False):
    if decode_only and (verify_only or recover_decoding or retry_memory_limits
                        or stop_after not in (None, 'decoding')):
        raise ValueError('Fresh decoding cannot select another pass or retry policy')
    if recover_decoding and (verify_only or retry_memory_limits or limit is not None or stop_after is not None):
        raise ValueError('Decoding recovery cannot select another pass, limit or retry policy')
    if verify_only and not resume:
        raise ValueError('Verification-only execution requires --resume')
    if verify_only and (retry_memory_limits or stop_after not in (None, 'verification')):
        raise ValueError('Verification-only execution cannot retry or select another pass')
    with ExitStack() as owners:
        output = Path(output)
        if recover_decoding:
            if (output/'decoding-recovery.json').exists():
                raise ValueError('Recovery output already has a registered successor; use that latest recovery')
            rows, requested = decoding_recovery_plan(index_path, output, memory_bytes, owners)
            registration = Path(index_path)/'decoding-recovery.json'
            registered = recovery_registration(output, requested)
            if registration.exists():
                if not resume or json.loads(registration.read_text()) != registered:
                    raise ValueError('Source already has a registered decoding recovery; resume its unchanged output')
            elif resume:
                raise ValueError('Decoding recovery resume requires its original registration')
            if resume:
                require_drained(output)
        else:
            raw = Path(index_path).read_bytes()
            index = json.loads(raw)
            if 'decoding_source' in index and not decode_only:
                raise ValueError('A saved catalogue index requires --decode-only')
            config = AAMSearchConfig(**index['config'])
            decoding_config = decoding_policy(index)
            verification = index.get('verification', 'enzymemap:verify')
            reference_policy = index.get('reference_policy',
                'exact annotated heavy-atom indices; no endpoint-symmetry relaxation')
            rows = [row for row in index['directions'] if limit is None or row['case'] < limit]
            requested = dict(input_index=str(Path(index_path).resolve()), input_index_sha256=hashlib.sha256(raw).hexdigest(),
                config=asdict(config), decode_config=asdict(decoding_config),
                cases=len({r['case'] for r in rows}), directions=len(rows),
                reference_policy=reference_policy)
            if decode_only:
                if not rows or not isinstance(index.get('decoding_source'), dict) or not index['decoding_source']:
                    raise ValueError('Fresh decoding requires saved catalogue provenance and nonempty scope')
                catalogues = [dict(case=row['case'], reversed=row['reversed'],
                                   catalogue=asdict(Blob(**row['catalogue']))) for row in rows]
                requested.update(mode='decode_only', decoding_source=index['decoding_source'],
                    decoder_source=code_identity(),
                    directions_to_decode=catalogues, memory_bytes=memory_bytes, case_seconds=case_seconds)
        identities = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}" for row in rows}
        if len(identities) != len(rows):
            raise ValueError('Duplicate reaction/direction inputs')
        requested = json.loads(json.dumps(requested))
        output.mkdir(parents=True, exist_ok=True)
        plan = output/'plan.json'
        if resume:
            saved_plan = json.loads(plan.read_text())
            if not recover_decoding:
                saved_plan['input_index'] = str(Path(saved_plan['input_index']).resolve())
            if saved_plan != requested:
                raise ValueError('Resume input/configuration differs')
            if decode_only:
                require_drained(output)
        else:
            if plan.exists():
                raise ValueError('Existing run requires --resume')
            plan.write_text(json.dumps(requested,indent=2)+'\n')
            if recover_decoding:
                pending = registration.with_suffix('.pending.json')
                with pending.open('w') as stream:
                    json.dump(registered, stream, indent=2)
                    stream.write('\n')
                    stream.flush()
                    os.fsync(stream.fileno())
                pending.replace(registration)
        stop = False
        def stop_feeding(signum, frame):
            nonlocal stop
            stop = True
        if not local:
            signal.signal(signal.SIGUSR1, stop_feeding)
        started = time.monotonic()
        status = dict(job=os.environ.get('SLURM_JOB_ID'), source=code_identity(), memory_bytes=memory_bytes,
                      expected_cases=requested['cases'], expected_directions=len(rows))
        if recover_decoding:
            status['decoding_recovery'] = True
        if decode_only:
            status['decoding_only'] = True
        if verify_only:
            if json.loads((output/'progress.json').read_text())['source'] != status['source']:
                raise ValueError('MAPPA source or compiled engine differs from the saved campaign')
            require_drained(output)
            status['verification_only'] = True
        stages = ([('decoding', None)] if recover_decoding or decode_only else
                  [('verification', 'deduplication')] if verify_only else
                  [('matching', None), ('deduplication', 'matching'),
                   ('decoding', 'deduplication'), ('verification', 'deduplication')])
        with tempfile.TemporaryDirectory(prefix='mappa-campaign-', dir=os.environ.get('SLURM_TMPDIR')) as scratch:
            ledger = Ledger(output, scratch, resume=resume)
            if verify_only:
                for name in identities:
                    prior = ledger.get(name, 'deduplication')
                    if (prior is None or prior[0] not in ('complete', 'error') or
                            (prior[0] == 'complete' and prior[1] is None)):
                        ledger.db.close()
                        raise ValueError('Verification requires terminal reduction for every direction '
                                         'and a reference for each completed catalogue')
            store = BundleStore(output/'bundles', scratch=scratch)
            factory = WorkPool.local if local else WorkPool.slurm
            with factory(workers=workers, store=store, memory_bytes=memory_bytes,
                         case_seconds=case_seconds, scratch=scratch) as pool:
                status['startup_seconds'] = pool.startup_elapsed_seconds
                status['workers'] = pool.queue.ready_workers()
                def checkpoint(stage, active):
                    ledger.checkpoint()
                    pool.checkpoint()
                    status.update(stage=stage, active_directions=len(active), completed=ledger.counts(),
                                  queue=pool.queue.snapshot(), elapsed_seconds=time.monotonic()-started)
                    pending = output/'progress.pending.json'
                    pending.write_text(json.dumps(status,indent=2)+'\n')
                    pending.replace(output/'progress.json')
                    print(json.dumps(status), flush=True)
                for stage, previous in stages:
                    active = {}; pending_rows = iter(rows); exhausted = False
                    checkpoint(stage, active)
                    next_checkpoint = time.monotonic()+30
                    while active or not exhausted:
                        stop |= time.monotonic()-started >= run_seconds
                        while not stop and not exhausted and len(active) < window:
                            try:
                                row = next(pending_rows)
                            except StopIteration:
                                exhausted = True
                                break
                            name = f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}"
                            existing = ledger.get(name, stage)
                            if existing:
                                if retry_memory_limits and stage == 'matching' and existing[0] == 'error':
                                    audit = ledger.retry_audit(name, stage)
                                    if audit is not None:
                                        record = store.load(audit)
                                        failed = [task for task in record['tasks'] if task[4] != 'complete']
                                        if not failed or any(task[3] >= memory_bytes for task in failed):
                                            continue
                                        try:
                                            result = pool.retry_matching(name, audit,
                                                states=('memory_limit', 'timeout'), retry_seconds=case_seconds)
                                        except ValueError:
                                            continue  # Other failure types still block the next pass.
                                        ledger.retire_failure(name, stage,
                                            f'resource retry at {memory_bytes} bytes with {case_seconds} fresh seconds')
                                        active[name] = result
                                continue
                            try:
                                if recover_decoding:
                                    result = pool.retry_decoding(name, Blob(**row['audit']),
                                                                 catalogue=Blob(**row['catalogue']))
                                elif decode_only:
                                    saved = StageResults('deduplication', store,
                                        references={name: Blob(**row['catalogue'])})
                                    result = pool.decode(saved, decoding_config)
                                elif stage == 'matching':
                                    result = pool.search({name: prepared_problem(row)}, config, seed_batch=1)
                                else:
                                    prior = ledger.get(name, previous)
                                    if prior is None or prior[0] != 'complete':
                                        raise RuntimeError('Previous pass did not complete')
                                    prior_results = StageResults(previous, store, references={name:Blob(**json.loads(prior[1]))})
                                    if stage == 'deduplication':
                                        result = pool.deduplicate(prior_results)
                                    elif stage == 'decoding':
                                        result = pool.decode(prior_results, decoding_config)
                                    else:
                                        catalogue = ledger.get(name, 'deduplication')
                                        key = f'{uuid.uuid4().hex}/{name}/verification'
                                        payload = dict(catalogue=Blob(**json.loads(catalogue[1])),
                                            reference=Blob(**row['reference']), reversed=row['reversed'])
                                        pool.queue.add([Task(key,verification,payload,case=key,
                                                            case_seconds=case_seconds,memory_bytes=memory_bytes)])
                                        result = StageResults(stage,store,pool=pool,tasks={name:(key,)})
                                active[name] = result
                            except Exception as error:
                                ledger.put(name,stage,'error',error=repr(error))
                        all_keys = [key for result in active.values() for group in result.tasks.values() for key in group]
                        states = pool.queue.states(all_keys)
                        for name,result in list(active.items()):
                            keys = result.tasks[name]
                            if any(states[key] in ('pending','running') for key in keys):
                                continue
                            failed = {key:states[key] for key in keys if states[key]!='complete'}
                            reference = None if failed else result.references[name]
                            # Publication precedes retirement; raw histories and attempts remain in bundles.
                            audit = pool.queue.archive_case(keys[0],store)
                            ledger.put(name,stage,'error' if failed else 'complete',reference,audit,failed)
                            del active[name]
                        if time.monotonic() >= next_checkpoint:
                            checkpoint(stage,active)
                            next_checkpoint = time.monotonic()+30
                        if stop and not active:
                            checkpoint(stage,active)
                            status['allocation_drained'] = True
                            (output/'progress.json').write_text(json.dumps(status,indent=2)+'\n')
                            return status
                        if active:
                            pool.check_workers()
                            time.sleep(.1)
                    failed = ledger.db.execute("SELECT count(*) FROM results WHERE stage=? AND state!='complete'",(stage,)).fetchone()[0]
                    if failed:
                        status.setdefault('failed_stages', []).append(stage)
                    checkpoint(stage,active)
                    if failed and (stage != 'decoding' or stop_after == 'decoding' or recover_decoding or decode_only):
                        raise RuntimeError(f'{stage} contains failed directions; next pass was not started')
                    if stage == 'verification' and not verify_only and status.get('failed_stages'):
                        raise RuntimeError('decoding contains failed directions; external verification completed')
                    if stage == stop_after:
                        status['pass_finished'] = stage
                        checkpoint(stage,{})
                        return status
                if verify_only or recover_decoding or decode_only:
                    stage = 'verification' if verify_only else 'decoding'
                    status['pass_finished'] = stage
                    checkpoint(stage,{})
                else:
                    status['finished'] = True
                    checkpoint('finished',{})
            ledger.db.close()
        return status


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('index'); parser.add_argument('output')
    parser.add_argument('--workers',type=int,default=32)
    parser.add_argument('--memory-gib',type=float,default=4.,help='Memory limit per compute worker (GiB)')
    parser.add_argument('--window',type=int,default=128)
    parser.add_argument('--limit',type=int)
    parser.add_argument('--local',action='store_true')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--verify-only',action='store_true',
                        help='Verify completed catalogues without requiring decoding; requires --resume')
    parser.add_argument('--recover-decoding',action='store_true',
                        help='Read an original campaign or its latest terminal recovery and retry decoding memory failures into separate output, preserving guards and remaining budgets')
    parser.add_argument('--decode-only',action='store_true',
                        help='Decode every saved catalogue in a prepared index into fresh output, with fresh recorded budgets')
    parser.add_argument('--retry-memory-limits',action='store_true',
                        help='Retry memory/time failures from smaller memory allocations with fresh time; preserve completed cuts and attempts')
    parser.add_argument('--stop-after',choices=['matching','deduplication','decoding','verification'])
    parser.add_argument('--run-seconds',type=float,default=12600.)
    parser.add_argument('--case-seconds',type=float,default=3600.,
                        help='Compute budget per fresh direction task; archived recovery retains its original remaining budget')
    args=parser.parse_args()
    status = run(args.index,args.output,workers=args.workers,memory_bytes=int(args.memory_gib*2**30),window=args.window,limit=args.limit,
        local=args.local,resume=args.resume,run_seconds=args.run_seconds,stop_after=args.stop_after,
        retry_memory_limits=args.retry_memory_limits,verify_only=args.verify_only,recover_decoding=args.recover_decoding,
        decode_only=args.decode_only,case_seconds=args.case_seconds)
    if status.get("allocation_drained"):
        raise SystemExit(75)


if __name__ == '__main__':
    main()
