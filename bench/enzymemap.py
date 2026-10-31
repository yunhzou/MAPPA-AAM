"""Fresh EnzymeMap passes over prepared inputs, using one persistent MAPPA pool."""
import argparse
from dataclasses import asdict
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

from mappa import AAMSearchConfig, query_path
from mappa.execution import BundleStore, StageResults, WorkPool
from mappa.execution.identity import code_identity
from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
from mappa.mappa_core.utils.io.records import problem_from_record


def prepared_problem(row):
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


def run(index_path, output, *, workers=32, window=128, memory_bytes=4*2**30,
        case_seconds=3600., local=False, limit=None, resume=False, run_seconds=12600., stop_after=None,
        retry_memory_limits=False):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    raw = Path(index_path).read_bytes()
    index = json.loads(raw)
    config = AAMSearchConfig(**index['config'])
    rows = [row for row in index['directions'] if limit is None or row['case'] < limit]
    identities = {f"{row['case']:05d}/{'reverse' if row['reversed'] else 'forward'}" for row in rows}
    if len(identities) != len(rows):
        raise ValueError('Duplicate reaction/direction inputs')
    requested = dict(input_index=str(Path(index_path).resolve()), input_index_sha256=hashlib.sha256(raw).hexdigest(),
        config=asdict(config), cases=len({r['case'] for r in rows}), directions=len(rows),
        reference_policy='exact annotated heavy-atom indices; no endpoint-symmetry relaxation')
    requested = json.loads(json.dumps(requested))
    plan = output/'plan.json'
    if resume:
        saved_plan = json.loads(plan.read_text())
        saved_plan['input_index'] = str(Path(saved_plan['input_index']).resolve())
        if saved_plan != requested:
            raise ValueError('Resume input/configuration differs')
    else:
        if plan.exists():
            raise ValueError('Existing run requires --resume')
        plan.write_text(json.dumps(requested,indent=2)+'\n')
    stop = False
    def stop_feeding(signum, frame):
        nonlocal stop
        stop = True
    if not local:
        signal.signal(signal.SIGUSR1, stop_feeding)
    started = time.monotonic()
    status = dict(job=os.environ.get('SLURM_JOB_ID'), source=code_identity(), memory_bytes=memory_bytes,
                  expected_cases=requested['cases'], expected_directions=len(rows))
    with tempfile.TemporaryDirectory(prefix='mappa-campaign-', dir=os.environ.get('SLURM_TMPDIR')) as scratch:
        ledger = Ledger(output, scratch, resume=resume)
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
            for stage, previous in [('matching',None),('deduplication','matching'),
                                    ('decoding','deduplication'),('verification','decoding')]:
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
                                    try:
                                        result = pool.retry_matching(name, audit)
                                    except ValueError:
                                        continue  # Other failure types still block the next pass.
                                    ledger.retire_failure(name, stage, 'explicit retry of diagnosed matching memory failure')
                                    active[name] = result
                            continue
                        try:
                            if stage == 'matching':
                                result = pool.search({name: prepared_problem(row)}, config, seed_batch=1)
                            else:
                                prior = ledger.get(name, previous)
                                if prior is None or prior[0] != 'complete':
                                    raise RuntimeError('Previous pass did not complete')
                                prior_results = StageResults(previous, store, references={name:Blob(**json.loads(prior[1]))})
                                if stage == 'deduplication':
                                    result = pool.deduplicate(prior_results)
                                elif stage == 'decoding':
                                    result = pool.decode(prior_results)
                                else:
                                    catalogue = ledger.get(name, 'deduplication')
                                    key = f'{uuid.uuid4().hex}/{name}/verification'
                                    payload = dict(catalogue=Blob(**json.loads(catalogue[1])),
                                        reference=Blob(**row['reference']), reversed=row['reversed'])
                                    pool.queue.add([Task(key,'enzymemap:verify',payload,case=key,
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
                checkpoint(stage,active)
                if ledger.db.execute("SELECT count(*) FROM results WHERE stage=? AND state!='complete'",(stage,)).fetchone()[0]:
                    raise RuntimeError(f'{stage} contains failed directions; next pass was not started')
                if stage == stop_after:
                    status['pass_finished'] = stage
                    checkpoint(stage,{})
                    return status
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
    parser.add_argument('--retry-memory-limits',action='store_true',
                        help='Retry archived matching memory failures; preserve completed cuts and attempts')
    parser.add_argument('--stop-after',choices=['matching','deduplication','decoding','verification'])
    parser.add_argument('--run-seconds',type=float,default=12600.)
    args=parser.parse_args()
    status = run(args.index,args.output,workers=args.workers,memory_bytes=int(args.memory_gib*2**30),window=args.window,limit=args.limit,
        local=args.local,resume=args.resume,run_seconds=args.run_seconds,stop_after=args.stop_after,
        retry_memory_limits=args.retry_memory_limits)
    if status.get("allocation_drained"):
        raise SystemExit(75)


if __name__ == '__main__':
    main()
