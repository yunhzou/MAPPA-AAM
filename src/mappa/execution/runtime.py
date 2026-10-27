"""The same supervised worker consumes an in-process or remote task queue."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import importlib
import multiprocessing as mp
from multiprocessing.managers import BaseManager
import os
from pathlib import Path
import signal
import tempfile
import threading
import time
import traceback
import uuid

from .bundles import BundleWriter, read_blob


class Workspace:
    """One prepared case per worker, plus bounded native caches inside that case."""
    def __init__(self, scratch, bundles):
        self.scratch = scratch
        self.writer = BundleWriter(bundles, scratch)
        self.prepared_key, self.prepared_value = None, None
        self.dependencies = {}

    def prepared(self, key, factory):
        if key != self.prepared_key:
            self.prepared_value = None
            self.prepared_key = None
            self.prepared_value = factory()
            self.prepared_key = key
        return self.prepared_value

    def load(self, blob):
        return read_blob(blob, self.scratch)

    def dependency(self, key):
        return self.load(self.dependencies[key]['blob'])


def _compute(connection, scratch, bundles):
    os.setsid()
    workspace = Workspace(scratch, bundles)
    connection.send(('ready', None))
    while True:
        work = connection.recv()
        if work is None:
            break
        task, dependencies = work
        workspace.dependencies = dependencies
        cpu, wall = time.process_time(), time.perf_counter()
        try:
            module, name = task.operation.split(':', 1)
            operation = getattr(importlib.import_module(module), name)
            output = operation(task.payload, workspace)
            metrics = dict(cpu_seconds=time.process_time()-cpu, compute_seconds=time.perf_counter()-wall)
            # Stop the work deadline at the compute boundary; allow evidence to publish.
            connection.send(('publishing', metrics))
            start = time.perf_counter()
            blob = workspace.writer.put(output)
            del output
            metrics['publish_seconds'] = time.perf_counter()-start
            connection.send(('complete', dict(blob=blob, **metrics)))
        except Exception:
            connection.send(('error', dict(error=traceback.format_exc(),
                cpu_seconds=time.process_time()-cpu, compute_seconds=time.perf_counter()-wall)))
        finally:
            workspace.dependencies = {}
    connection.close()


def _stop(process):
    if process.is_alive():
        # Child operations may not start pools. Kill their process group if a
        # library nevertheless spawned children, preserving the stated guard.
        try:
            os.killpg(process.pid, signal.SIGTERM) if os.getpgid(process.pid) == process.pid else process.terminate()
        except ProcessLookupError:
            pass
        process.join(2)
        if process.is_alive():
            os.killpg(process.pid, signal.SIGKILL)
    process.join()


def worker(queue, *, identity, scratch, bundles, memory_bytes, poll_seconds=.2, publish_seconds=120.):
    """One compute CPU, no nested pool; limits produce explicit incomplete states.

    psutil is required by the execution extra. The task's memory request is an
    admission estimate; memory_bytes is this worker's actual process-tree guard.
    """
    import psutil
    worker_id = f'{os.uname().nodename}:{os.getpid()}:{uuid.uuid4().hex}'
    with tempfile.TemporaryDirectory(prefix='mappa-worker-', dir=scratch) as temporary:
        parent, child = mp.get_context('spawn').Pipe()
        process = mp.get_context('spawn').Process(target=_compute, args=(child, temporary, bundles))
        process.start()
        child.close()
        try:
            if not parent.poll(60) or parent.recv()[0] != 'ready':
                raise RuntimeError('Compute worker failed to initialize')
            while True:
                claim = queue.claim(worker_id, memory_bytes, identity)
                if claim['state'] in ('finished', 'insufficient_memory'):
                    return claim['state']
                if claim['state'] == 'waiting':
                    time.sleep(max(.5, poll_seconds))
                    continue
                task, token = claim['task'], claim['token']
                parent.send((task, claim['dependencies']))
                peak, publishing, publication_deadline = 0, False, None
                start, heartbeat = time.perf_counter(), time.monotonic()
                result = {}
                while True:
                    try:
                        tree = psutil.Process(process.pid)
                        rss = sum(p.memory_info().rss for p in [tree, *tree.children(recursive=True)])
                        peak = max(peak, rss)
                    except psutil.NoSuchProcess:
                        continue
                    if parent.poll(poll_seconds):
                        state, result = parent.recv()
                        if state == 'publishing':
                            if time.time() > claim['deadline']:
                                state = 'timeout'
                                break
                            publishing = True
                            publication_deadline = time.monotonic()+publish_seconds
                            continue
                        break
                    if not process.is_alive():
                        state, result = 'error', dict(exit_code=process.exitcode)
                        break
                    if peak > memory_bytes:
                        state = 'memory_limit'
                        break
                    if not publishing and time.time() >= claim['deadline']:
                        state = 'timeout'
                        break
                    if publishing and time.monotonic() >= publication_deadline:
                        state, result = 'error', dict(error='Evidence publication deadline exceeded')
                        break
                    if time.monotonic()-heartbeat >= 5:
                        if not queue.heartbeat(task.key, token):
                            return  # lease lost; never accept stale work
                        heartbeat = time.monotonic()
                result.update(peak_tree_rss_bytes=peak, worker=worker_id,
                              elapsed_seconds=time.perf_counter()-start)
                queue.finish(task.key, token, state, result)
                if state != 'complete':
                    # Reset only the worker after a failed computation, never
                    # resubmit the task or substitute a different algorithm.
                    return 'failed'
        except (EOFError, BrokenPipeError):
            if 'task' in locals():
                queue.finish(task.key, token, 'error', dict(error='Compute process exited', exit_code=process.exitcode))
                return 'failed'
            raise
        finally:
            _stop(process)
            parent.close()


def _supervise(queue, options):
    while worker(queue, **options) == 'failed':
        pass


def run(queue, *, workers=1, identity, scratch, bundles, memory_bytes):
    """Use the common worker with one CPU or several local compute processes.

    memory_bytes is a per-worker reservation: workers * memory_bytes must fit
    the allocation. Slurm workers use this same entry point against connect().
    """
    if workers < 1 or memory_bytes < 1:
        raise ValueError('Positive worker count and per-worker memory are required')
    options = dict(identity=identity, scratch=scratch, bundles=bundles, memory_bytes=memory_bytes)
    if workers == 1:
        _supervise(queue, options)
    else:
        with ThreadPoolExecutor(max_workers=workers) as supervisors:
            futures = [supervisors.submit(_supervise, queue, options) for _ in range(workers)]
            for future in futures:
                future.result()
    return queue.snapshot()


@contextmanager
def coordinator(queue, address, authkey):
    """Authenticated trusted-cluster transport; metadata stays with its owner.

    Use a private network or SSH tunnel. Never expose Python manager/pickle
    endpoints to untrusted clients. Credentials belong outside source control.
    """
    class Manager(BaseManager):
        pass
    Manager.register('queue', callable=lambda: queue,
                     exposed=('claim', 'heartbeat', 'finish', 'snapshot', 'results'))
    server = Manager(address=address, authkey=authkey).get_server()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server.address
    finally:
        server.stop_event.set()
        server.listener.close()
        thread.join(2)


def connect(address, authkey):
    class Manager(BaseManager):
        pass
    Manager.register('queue')
    manager = Manager(address=address, authkey=authkey)
    manager.connect()
    return manager.queue()
