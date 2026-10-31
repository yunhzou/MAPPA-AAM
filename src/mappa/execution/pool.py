"""Persistent execution around separate scientific stages; no matching logic."""
from contextlib import ExitStack
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import time
import uuid

from mappa.mappa_core.core.scheduler import Task, TaskQueue
from mappa.mappa_core.core.scheduler.runtime import run, coordinator, connect_store
from mappa.mappa_core.core.scheduler.search import plan_search
from mappa.mappa_core.core.scheduler.storage import MemoryStore, BundleStore
from mappa.mappa_core.core.scheduler.bundles import Blob
from mappa.mappa_core.core.scheduler.results import MatchingResult, CutScope
from mappa.mappa_core.domain import AAMSearchConfig
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
from .identity import code_identity


class StageResults:
    """Named typed artifacts; loading a reaction is explicit and lazy."""
    def __init__(self, stage, store, *, pool=None, tasks=None, plans=None, references=None):
        self.stage, self.store, self.pool = stage, store, pool
        self.tasks = tasks or {}
        self.plans = plans or {}
        self._references = references

    def _rows(self):
        keys = [key for group in self.tasks.values() for key in group]
        return self.pool.queue.results(keys)

    def wait_and_validate(self):
        if self._references is not None:
            return self
        while True:
            rows = self._rows()
            failed = {key: row['state'] for key, row in rows.items()
                      if row['state'] not in ('pending', 'running', 'complete')}
            if failed:
                raise RuntimeError(f'{self.stage} did not complete: {failed}')
            if all(row['state'] == 'complete' for row in rows.values()):
                self._collect(rows)
                self.pool.checkpoint()
                return self
            self.pool.check_workers()
            time.sleep(.01)

    def _collect(self, rows):
        self._references = {}
        for case, keys in self.tasks.items():
            references = tuple(rows[key]['result']['blob'] for key in keys)
            if self.stage == 'matching':
                source, scope = self.plans[case]
                reference = self.store.put(MatchingResult(case, source, scope, references))
            else:
                reference, = references
            self._references[case] = reference

    @property
    def references(self):
        if self._references is None:
            rows = self._rows()
            if any(row['state'] != 'complete' for row in rows.values()):
                raise RuntimeError('Finish and validate the previous stage before submitting the next')
            self._collect(rows)
        return dict(self._references)

    def load(self, case):
        return self.store.load(self.references[case])

    def __iter__(self):
        return iter(self.tasks if self._references is None else self._references)

    def save_manifest(self, path):
        references = self.references
        if not all(isinstance(ref, Blob) for ref in references.values()):
            raise ValueError('A durable manifest requires BundleStore; RAM references are process-lifetime only')
        value = dict(schema='mappa.stage/v1', stage=self.stage,
                     source=code_identity(),
                     results={case: asdict(ref) for case, ref in references.items()})
        path = Path(path)
        temporary = path.with_suffix(path.suffix+'.pending')
        temporary.write_text(json.dumps(value, indent=2)+'\n')
        temporary.replace(path)

    @classmethod
    def open(cls, path, *, store):
        value = json.loads(Path(path).read_text())
        if value['schema'] != 'mappa.stage/v1' or value['stage'] not in ('matching', 'deduplication', 'decoding'):
            raise ValueError('Unsupported stage manifest')
        return cls(value['stage'], store,
                   references={case: Blob(**ref) for case, ref in value['results'].items()})


class WorkPool:
    """One worker lifetime, explicit stages, RAM or aggregated bundle transport."""
    def __init__(self, *, workers=1, store=None, memory_bytes=2**30,
                 case_seconds=3600., scratch=None, slurm=False, startup_seconds=120.):
        if workers < 1 or memory_bytes < 1 or case_seconds <= 0 or startup_seconds <= 0:
            raise ValueError('Positive workers, memory and stage time are required')
        self.workers = workers
        self.store = store if store is not None else MemoryStore()
        self.memory_bytes, self.case_seconds = memory_bytes, case_seconds
        self.startup_seconds = startup_seconds
        self.scratch_parent, self.use_slurm = scratch, slurm
        self.identity = uuid.uuid4().hex
        self.code_identity = code_identity()
        self._stack = None
        self._runner = self._launch = None

    @classmethod
    def local(cls, workers=1, **options):
        return cls(workers=workers, **options)

    @classmethod
    def slurm(cls, *, workers=None, **options):
        if 'SLURM_JOB_ID' not in os.environ:
            raise RuntimeError('WorkPool.slurm requires an existing Slurm allocation')
        cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
        workers = cpus if workers is None else workers
        if workers > cpus:
            raise ValueError('Workers exceed CPUs allocated to each node task')
        return cls(workers=workers, slurm=True, **options)

    def __enter__(self):
        started = time.monotonic()
        self._stack = ExitStack()
        try:
            self.scratch = Path(self._stack.enter_context(
                tempfile.TemporaryDirectory(prefix='mappa-pool-', dir=self.scratch_parent)))
            database = self.scratch/'queue.sqlite' if isinstance(self.store, BundleStore) else ':memory:'
            self.queue = TaskQueue(database, self.identity, keep_open=True)
            self._stack.callback(self.queue.close)
            worker_store = self.store
            if self.use_slurm or isinstance(self.store, MemoryStore):
                key = os.urandom(32)
                host = socket.getfqdn() if self.use_slurm else '127.0.0.1'
                address = self._stack.enter_context(coordinator(self.queue, (host, 0), key,
                    store=self.store if isinstance(self.store, MemoryStore) else None))
                if isinstance(self.store, MemoryStore):
                    worker_store = connect_store(address, key)
            if self.use_slurm:
                self._start_slurm(address, key)
            else:
                supervisors = self._stack.enter_context(ThreadPoolExecutor(max_workers=1))
                self._runner = supervisors.submit(run, self.queue, workers=self.workers,
                    identity=self.identity, scratch=self.scratch, store=worker_store,
                    memory_bytes=self.memory_bytes)
            expected = self.workers * (int(os.environ['SLURM_JOB_NUM_NODES']) if self.use_slurm else 1)
            while len(self.queue.ready_workers()) < expected:
                self.check_workers()
                if time.monotonic()-started >= self.startup_seconds:
                    raise TimeoutError(f'Only {len(self.queue.ready_workers())}/{expected} workers became ready')
                time.sleep(.01)
            self.startup_elapsed_seconds = time.monotonic()-started
            return self
        except BaseException:
            if hasattr(self, 'queue'):
                self.queue.finish_input(cancel_pending=True)
            if self._launch is not None and self._launch.poll() is None:
                self._launch.terminate()
                self._launch.wait()
            self._stack.close()
            raise

    def _start_slurm(self, address, key):
        nodes = int(os.environ['SLURM_JOB_NUM_NODES'])
        cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
        allocated_memory = os.environ.get('SLURM_MEM_PER_NODE')
        if allocated_memory and self.workers*self.memory_bytes > int(allocated_memory)*2**20:
            raise ValueError('Worker memory reservations exceed the Slurm node allocation')
        command = ['srun', '--nodes', str(nodes), '--ntasks', str(nodes),
                   '--ntasks-per-node', '1', '--cpus-per-task', str(cpus),
                   sys.executable, '-m', 'mappa.execution',
                   '--host', address[0], '--port', str(address[1]),
                   '--identity', self.identity, '--workers', str(self.workers),
                   '--code-identity', self.code_identity,
                   '--memory-bytes', str(self.memory_bytes)]
        if isinstance(self.store, BundleStore):
            command.extend(('--bundles', str(self.store.directory)))
        self._launch = subprocess.Popen(command, env={**os.environ, 'MAPPA_POOL_AUTH': key.hex()})

    def check_workers(self):
        if self._runner is not None and self._runner.done():
            self._runner.result()
            raise RuntimeError('Workers stopped before the requested stage completed')
        if self._launch is not None and self._launch.poll() is not None:
            raise RuntimeError(f'Slurm workers exited with code {self._launch.returncode}')

    def checkpoint(self):
        if isinstance(self.store, BundleStore):
            self.queue.backup(self.store.directory/'queue.sqlite')

    def search(self, problems, config=None, *, seed_batch=1, execution='reused_native'):
        """Submit named reaction/direction inputs; only matching runs here."""
        config = config or AAMSearchConfig()
        keys, plans, tasks = {}, {}, []
        batch = uuid.uuid4().hex
        for case, problem in problems.items():
            if not isinstance(case, str) or not case:
                raise ValueError('Cases require nonempty string identities')
            planned, cuts = plan_search(f'{batch}/{case}', problem, config, self.store,
                seed_batch=seed_batch, case_seconds=self.case_seconds,
                memory_bytes=self.memory_bytes, execution=execution)
            finals = {task.key: task for task in planned}
            scopes = tuple(CutScope(finals[key].payload['number'], tuple(finals[key].payload['cut']),
                tuple(range(finals[key].payload['seeds']))) for key in cuts)
            plans[case] = planned[0].payload['input'], scopes
            keys[case] = cuts
            tasks.extend(planned)
        self.queue.add(tasks)
        self.queue.optimize()
        return StageResults('matching', self.store, pool=self, tasks=keys, plans=plans)

    def retry_matching(self, case, audit, *, states=('memory_limit',)):
        """Retry diagnosed failures; reuse completed cuts and retain every attempt."""
        record = self.store.load(audit)
        failed = [row for row in record['tasks'] if row[4] != 'complete']
        if not failed or any(row[4] not in (*states, 'blocked') for row in failed):
            raise ValueError('Audit contains outcomes outside the requested retry states')
        tasks = [_loads_checkpoint(row[1]) for row in record['tasks']]
        prefix = 'mappa.mappa_core.core.scheduler.search:'
        if any(task.operation not in (prefix+'search_batch', prefix+'finalize_cut', prefix+'search_complete_cut') for task in tasks):
            raise ValueError('Matching recovery requires a matching task audit')
        finals = sorted((task for task in tasks if task.operation != prefix+'search_batch'),
                        key=lambda task: task.payload['number'])
        scopes = tuple(CutScope(task.payload['number'], tuple(task.payload['cut']),
                               tuple(range(task.payload['seeds']))) for task in finals)
        self.queue.restore_case(record)
        for row in failed:
            if row[4] in states:
                self.queue.retry(row[0])
        self.queue.optimize()
        return StageResults('matching', self.store, pool=self,
            tasks={case: tuple(task.key for task in finals)},
            plans={case: (tasks[0].payload['input'], scopes)})

    def _stage(self, previous, expected, stage, operation, field, **options):
        if previous.stage != expected:
            raise TypeError(f'{stage} requires {expected} results')
        if type(previous.store) is not type(self.store):
            raise ValueError('Stages must use the same storage mode')
        if isinstance(previous.store, MemoryStore) and previous.store is not self.store:
            raise ValueError('RAM stages must use the same MemoryStore')
        references = previous.references  # Global pass barrier, including failed/missing work.
        tasks, keys = [], {}
        batch = uuid.uuid4().hex
        for case, reference in references.items():
            key = f'{batch}/{case}/{stage}'
            tasks.append(Task(key, operation, {field: reference, **options},
                case=key, case_seconds=self.case_seconds, memory_bytes=self.memory_bytes))
            keys[case] = (key,)
        self.queue.add(tasks)
        self.queue.optimize()
        return StageResults(stage, self.store, pool=self, tasks=keys)

    def deduplicate(self, matched):
        return self._stage(matched, 'matching', 'deduplication',
            'mappa.mappa_core.core.scheduler.decoding:reduce_matching', 'matching')

    def decode(self, catalogues, config=None):
        return self._stage(catalogues, 'deduplication', 'decoding',
            'mappa.mappa_core.core.scheduler.decoding:decode_catalogue', 'catalogue',
            config=config or EventDecodeConfig())

    def __exit__(self, exc_type, exc, traceback):
        try:
            self.queue.finish_input(cancel_pending=exc is not None)
            if self._runner is not None:
                self._runner.result()
            if self._launch is not None:
                code = self._launch.wait()
                if code and exc is None:
                    raise RuntimeError(f'Slurm workers exited with code {code}')
            self.checkpoint()
        finally:
            self._stack.close()
