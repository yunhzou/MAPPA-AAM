"""One coordinator owns task metadata; workers never open its SQLite file."""
from dataclasses import dataclass
from graphlib import TopologicalSorter
import pickle
import os
import tempfile
from pathlib import Path
import sqlite3
import threading
import time
import uuid


@dataclass(frozen=True)
class Task:
    key: str
    operation: str
    payload: object
    dependencies: tuple[str, ...] = ()
    case: str = ''
    case_seconds: float = 3600.
    memory_bytes: int = 2**30


class TaskQueue:
    """Durable dependency queue for trusted workers running the same code.

    Place the database on the coordinator's local disk, with explicit backups
    to shared storage. Lost leases require an explicit retry; no silent reruns.
    """
    def __init__(self, path, identity, *, lease_seconds=120.):
        self.lock = threading.RLock()
        self.lease_seconds = lease_seconds
        self.db = sqlite3.connect(str(path), check_same_thread=False)
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.executescript('''
            CREATE TABLE IF NOT EXISTS identity (value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS cases (name TEXT PRIMARY KEY, seconds REAL, deadline REAL);
            CREATE TABLE IF NOT EXISTS tasks (
                key TEXT PRIMARY KEY, specification BLOB, case_name TEXT, memory INTEGER,
                state TEXT, worker TEXT, token TEXT, expires REAL, result BLOB);
            CREATE TABLE IF NOT EXISTS attempts (
                token TEXT PRIMARY KEY, task TEXT, worker TEXT, started REAL, finished REAL, state TEXT, result BLOB);
            CREATE TABLE IF NOT EXISTS dependencies (task TEXT, dependency TEXT, PRIMARY KEY(task, dependency));
            CREATE INDEX IF NOT EXISTS case_deadline ON cases(deadline);
            CREATE INDEX IF NOT EXISTS task_state ON tasks(state);
            CREATE INDEX IF NOT EXISTS task_case ON tasks(case_name,state);
            CREATE INDEX IF NOT EXISTS task_expiry ON tasks(state,expires);
            CREATE INDEX IF NOT EXISTS dependency_parent ON dependencies(dependency);
        ''')
        row = self.db.execute('SELECT value FROM identity').fetchone()
        if row and row[0] != identity:
            self.db.close()
            raise ValueError('Queue code/input identity differs')
        if row is None:
            self.db.execute('INSERT INTO identity VALUES (?)', (identity,))
            self.db.commit()
        self.identity = identity

    def add(self, tasks):
        tasks = list(tasks)
        with self.lock, self.db:
            if len({t.key for t in tasks}) != len(tasks):
                raise ValueError('Duplicate task key')
            new = {}
            for task in tasks:
                if task.case_seconds <= 0 or task.memory_bytes <= 0:
                    raise ValueError('Resource budgets must be positive')
                previous = self.db.execute('SELECT specification FROM tasks WHERE key=?', (task.key,)).fetchone()
                if previous:
                    if previous[0] != pickle.dumps(task, protocol=5):
                        raise ValueError('Task identity reused with different work')
                else:
                    new[task.key] = task
            external = {d for t in new.values() for d in t.dependencies if d not in new}
            for dependency in external:
                if self.db.execute('SELECT 1 FROM tasks WHERE key=?', (dependency,)).fetchone() is None:
                    raise ValueError('Unknown dependency')
            # Existing specifications are immutable and cannot refer to future
            # tasks. Only the new subgraph can introduce a cycle.
            tuple(TopologicalSorter({k: [d for d in t.dependencies if d in new]
                                     for k, t in new.items()}).static_order())
            for task in tasks:
                data = pickle.dumps(task, protocol=5)
                previous = self.db.execute('SELECT specification FROM tasks WHERE key=?', (task.key,)).fetchone()
                if previous:
                    if previous[0] != data:
                        raise ValueError('Task identity reused with different work')
                    continue
                case = task.case or task.key
                row = self.db.execute('SELECT seconds FROM cases WHERE name=?', (case,)).fetchone()
                if row and row[0] != task.case_seconds:
                    raise ValueError('Tasks in one case must share the same budget')
                self.db.execute('INSERT OR IGNORE INTO cases VALUES (?,?,NULL)', (case, task.case_seconds))
                self.db.execute('INSERT INTO tasks VALUES (?,?,?,?,?,?,?,?,?)',
                    (task.key, data, case, task.memory_bytes, 'pending', None, None, None, None))
                self.db.executemany('INSERT INTO dependencies VALUES (?,?)',
                    [(task.key, d) for d in task.dependencies])

    def _refresh(self, now):
        self.db.execute("""UPDATE attempts SET state='lost',finished=? WHERE state='running' AND token IN
            (SELECT token FROM tasks WHERE state='running' AND expires<?)""", (now, now))
        self.db.execute("UPDATE tasks SET state='lost' WHERE state='running' AND expires<?", (now,))
        self.db.execute("""UPDATE tasks SET state='timeout' WHERE state='pending' AND case_name IN
            (SELECT name FROM cases WHERE deadline IS NOT NULL AND deadline<=?)""", (now,))
        while self.db.execute("""UPDATE tasks SET state='blocked' WHERE state='pending' AND key IN
            (SELECT d.task FROM dependencies d JOIN tasks p ON p.key=d.dependency
             WHERE p.state IN ('error','timeout','memory_limit','lost','blocked'))""").rowcount:
            pass

    def claim(self, worker, memory_bytes, identity):
        if identity != self.identity:
            raise ValueError('Worker code/input identity differs')
        with self.lock, self.db:
            now = time.time()
            self._refresh(now)
            row = self.db.execute("""SELECT key,specification,case_name FROM tasks t
                WHERE state='pending' AND memory<=? AND NOT EXISTS (
                  SELECT 1 FROM dependencies d JOIN tasks p ON p.key=d.dependency
                  WHERE d.task=t.key AND p.state!='complete') ORDER BY t.rowid LIMIT 1""",
                (memory_bytes,)).fetchone()
            if row is None:
                active = self.db.execute("SELECT count(*) FROM tasks WHERE state IN ('pending','running')").fetchone()[0]
                if not active:
                    return {'state': 'finished'}
                running = self.db.execute("SELECT count(*) FROM tasks WHERE state='running'").fetchone()[0]
                return {'state': 'waiting' if running else 'insufficient_memory'}
            key, specification, case = row
            self.db.execute('UPDATE cases SET deadline=?+seconds WHERE name=? AND deadline IS NULL', (now, case))
            deadline = self.db.execute('SELECT deadline FROM cases WHERE name=?', (case,)).fetchone()[0]
            token = uuid.uuid4().hex
            self.db.execute("UPDATE tasks SET state='running',worker=?,token=?,expires=? WHERE key=?",
                            (worker, token, now+self.lease_seconds, key))
            self.db.execute('INSERT INTO attempts VALUES (?,?,?,?,NULL,?,NULL)',
                            (token, key, worker, now, 'running'))
            dependencies = {d: pickle.loads(r) for d, r in self.db.execute('''
                SELECT p.key,p.result FROM dependencies d JOIN tasks p ON p.key=d.dependency WHERE d.task=?''', (key,))}
            return dict(state='work', task=pickle.loads(specification), token=token,
                        deadline=deadline, dependencies=dependencies)

    def heartbeat(self, key, token):
        with self.lock, self.db:
            return bool(self.db.execute("UPDATE tasks SET expires=? WHERE key=? AND token=? AND state='running'",
                (time.time()+self.lease_seconds, key, token)).rowcount)

    def finish(self, key, token, state, result):
        if state not in ('complete', 'error', 'timeout', 'memory_limit', 'lost'):
            raise ValueError('Invalid task outcome')
        with self.lock, self.db:
            self._refresh(time.time())
            data = pickle.dumps(result, protocol=5)
            accepted = bool(self.db.execute("UPDATE tasks SET state=?,result=? WHERE key=? AND token=? AND state='running'",
                (state, data, key, token)).rowcount)
            if accepted:
                self.db.execute('UPDATE attempts SET state=?,result=?,finished=? WHERE token=?',
                                (state, data, time.time(), token))
            return accepted


    def retry(self, key):
        """Explicit recovery after diagnosis; completed work and budgets survive."""
        with self.lock, self.db:
            row = self.db.execute('SELECT state FROM tasks WHERE key=?', (key,)).fetchone()
            if row is None or row[0] not in ('lost', 'error', 'memory_limit'):
                raise ValueError('Only diagnosed failed/lost tasks may be retried')
            self.db.execute("UPDATE tasks SET state='pending',token=NULL,worker=NULL WHERE key=?", (key,))
            self.db.execute("UPDATE tasks SET state='pending' WHERE state='blocked'")
            self._refresh(time.time())

    def snapshot(self):
        with self.lock, self.db:
            self._refresh(time.time())
            return dict(self.db.execute('SELECT state,count(*) FROM tasks GROUP BY state'))

    def results(self, keys=None):
        with self.lock:
            rows = (self.db.execute('SELECT key,state,result FROM tasks ORDER BY rowid') if keys is None
                    else (self.db.execute('SELECT key,state,result FROM tasks WHERE key=?', (key,)).fetchone()
                          for key in keys))
            return {k: dict(state=s, result=pickle.loads(r) if r else None) for k, s, r in rows}

    def optimize(self):
        """Prepare planner statistics after bulk task insertion, before workers.

        Without cardinalities SQLite can scan every pending task for an empty
        dependency update. Analysis changes query plans, never queue outcomes.
        """
        with self.lock, self.db:
            self.db.execute('PRAGMA analysis_limit=1000')
            self.db.execute('ANALYZE')

    def backup(self, path):
        """One compact consistent checkpoint, never a shared database writer."""
        path = Path(path)
        with self.lock, tempfile.TemporaryDirectory(dir=path.parent) as temporary:
            snapshot = Path(temporary)/'queue.sqlite'
            destination = sqlite3.connect(snapshot)
            try:
                self.db.backup(destination)
                # The backup inherits WAL mode. Consolidate it before publishing
                # so copying this single file preserves every committed row.
                destination.execute('PRAGMA journal_mode=DELETE')
            finally:
                destination.close()
            os.replace(snapshot, path)

    def close(self):
        self.db.close()
