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
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint


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
        self.path = Path(path)
        self.db = sqlite3.connect(str(path), check_same_thread=False, timeout=60)
        columns = {r[1] for r in self.db.execute('PRAGMA table_info(cases)')}
        if columns and 'spent' not in columns:
            self.db.close()
            raise ValueError('Unsupported queue schema')
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.executescript('''
            CREATE TABLE IF NOT EXISTS identity (value TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS cases (
                name TEXT PRIMARY KEY, seconds REAL, deadline REAL,
                spent REAL NOT NULL DEFAULT 0, active_since REAL);
            CREATE TABLE IF NOT EXISTS active_work (token TEXT PRIMARY KEY, case_name TEXT);
            CREATE INDEX IF NOT EXISTS active_case ON active_work(case_name);
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
            CREATE TABLE IF NOT EXISTS resource_changes (
                utc REAL, reason TEXT, memory INTEGER, changed_tasks INTEGER);
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
            specifications = {}
            for task in tasks:
                if task.case_seconds <= 0 or task.memory_bytes <= 0:
                    raise ValueError('Resource budgets must be positive')
                data = pickle.dumps(task, protocol=5)
                previous = self.db.execute('SELECT specification FROM tasks WHERE key=?', (task.key,)).fetchone()
                if previous:
                    if previous[0] != data:
                        raise ValueError('Task identity reused with different work')
                else:
                    new[task.key] = task
                    specifications[task.key] = data
            external = {d for t in new.values() for d in t.dependencies if d not in new}
            for dependency in external:
                if self.db.execute('SELECT 1 FROM tasks WHERE key=?', (dependency,)).fetchone() is None:
                    raise ValueError('Unknown dependency')
            # Existing specifications are immutable and cannot refer to future
            # tasks. Only the new subgraph can introduce a cycle.
            tuple(TopologicalSorter({k: [d for d in t.dependencies if d in new]
                                     for k, t in new.items()}).static_order())
            for task in new.values():
                case = task.case or task.key
                row = self.db.execute('SELECT seconds FROM cases WHERE name=?', (case,)).fetchone()
                if row and row[0] != task.case_seconds:
                    raise ValueError('Tasks in one case must share the same budget')
                self.db.execute('INSERT OR IGNORE INTO cases(name,seconds) VALUES (?,?)', (case, task.case_seconds))
                self.db.execute('INSERT INTO tasks VALUES (?,?,?,?,?,?,?,?,?)',
                    (task.key, specifications[task.key], case, task.memory_bytes, 'pending', None, None, None, None))
                self.db.executemany('INSERT INTO dependencies VALUES (?,?)',
                    [(task.key, d) for d in task.dependencies])

    def _refresh(self, now):
        for token, expires in self.db.execute(
                "SELECT token,expires FROM tasks WHERE state='running' AND expires<? ORDER BY expires", (now,)).fetchall():
            self._end_work(token, expires)
        self.db.execute("""UPDATE attempts SET state='lost',finished=? WHERE state='running' AND token IN
            (SELECT token FROM tasks WHERE state='running' AND expires<?)""", (now, now))
        self.db.execute("UPDATE tasks SET state='lost' WHERE state='running' AND expires<?", (now,))
        self.db.execute("""UPDATE tasks SET state='timeout' WHERE state='pending' AND case_name IN
            (SELECT name FROM cases WHERE spent>=seconds OR deadline<=?)""", (now,))
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
            self.db.execute('UPDATE cases SET active_since=?,deadline=?+seconds-spent WHERE name=? AND active_since IS NULL',
                            (now, now, case))
            deadline = self.db.execute('SELECT deadline FROM cases WHERE name=?', (case,)).fetchone()[0]
            token = uuid.uuid4().hex
            self.db.execute('INSERT INTO active_work VALUES (?,?)', (token, case))
            self.db.execute("UPDATE tasks SET state='running',worker=?,token=?,expires=? WHERE key=?",
                            (worker, token, now+self.lease_seconds, key))
            self.db.execute('INSERT INTO attempts VALUES (?,?,?,?,NULL,?,NULL)',
                            (token, key, worker, now, 'running'))
            dependencies = {d: _loads_checkpoint(r) for d, r in self.db.execute('''
                SELECT p.key,p.result FROM dependencies d JOIN tasks p ON p.key=d.dependency WHERE d.task=?''', (key,))}
            return dict(state='work', task=_loads_checkpoint(specification), token=token,
                        deadline=deadline, dependencies=dependencies)

    def _end_work(self, token, now):
        row = self.db.execute('SELECT case_name FROM active_work WHERE token=?', (token,)).fetchone()
        if row is None:
            return
        case = row[0]
        self.db.execute('DELETE FROM active_work WHERE token=?', (token,))
        if self.db.execute('SELECT 1 FROM active_work WHERE case_name=? LIMIT 1', (case,)).fetchone() is None:
            self.db.execute('UPDATE cases SET spent=spent+max(0,?-active_since),active_since=NULL,deadline=NULL WHERE name=?',
                            (now, case))

    def computed(self, key, token):
        """Pause a direction once its last computation reaches publication.

        Overlapping tasks share one active wall clock. Queued dependencies,
        gaps between allocations and output publication do not consume it.
        """
        with self.lock, self.db:
            now = time.time()
            self._refresh(now)
            valid = self.db.execute("SELECT 1 FROM tasks WHERE key=? AND token=? AND state='running'",
                                    (key, token)).fetchone() is not None
            if valid:
                self._end_work(token, now)
            return valid

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
                self._end_work(token, time.time())
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

    def set_pending_memory(self, memory_bytes, *, reason):
        """Change admission estimates without changing saved work or remaining budgets.

        A node-local maintenance connection may call this during a worker
        handoff. SQLite serializes it with claims and publications. Setting an
        estimate above the old workers' capacity lets their active work finish
        before they retire through the normal insufficient-memory response.
        """
        if memory_bytes <= 0 or not reason.strip():
            raise ValueError('A positive estimate and recorded reason are required')
        with self.lock, self.db:
            self.db.execute('INSERT INTO resource_changes VALUES (?,?,?,?)',
                            (time.time(), reason, memory_bytes, 0))
            change = self.db.execute('SELECT last_insert_rowid()').fetchone()[0]
        changed, cursor = 0, 0
        while True:
            # Bound write-lock duration so a live coordinator can keep accepting
            # heartbeats and completed outputs during a large queue handoff.
            with self.lock, self.db:
                self.db.execute('BEGIN IMMEDIATE')
                rows = self.db.execute("SELECT rowid FROM tasks WHERE rowid>? AND state='pending' AND memory!=? ORDER BY rowid LIMIT 250",
                                       (cursor, memory_bytes)).fetchall()
                if not rows:
                    break
                cursor = rows[-1][0]
                changed += self.db.executemany("UPDATE tasks SET memory=? WHERE rowid=? AND state='pending'",
                                               [(memory_bytes, row[0]) for row in rows]).rowcount
                self.db.execute('UPDATE resource_changes SET changed_tasks=? WHERE rowid=?', (changed, change))
            time.sleep(.005)
        return changed

    def retry_memory_limits(self, previous_bytes, memory_bytes, *, reason):
        """Explicit larger-memory stage; preserve attempts and case deadlines."""
        if not 0 < previous_bytes < memory_bytes or not reason.strip():
            raise ValueError('The next stage must increase memory and record a reason')
        with self.lock, self.db:
            # Acquire the write transaction before selecting retry candidates.
            self.db.execute('INSERT INTO resource_changes VALUES (?,?,?,0)',
                            (time.time(), reason, memory_bytes))
            change = self.db.execute('SELECT last_insert_rowid()').fetchone()[0]
            keys = [key for key, raw in self.db.execute(
                "SELECT key,result FROM tasks WHERE state='memory_limit'")
                if _loads_checkpoint(raw).get('memory_guard_bytes') == previous_bytes]
            self.db.executemany("UPDATE tasks SET state='pending',memory=?,token=NULL,worker=NULL WHERE key=?",
                                [(memory_bytes, key) for key in keys])
            if keys:
                self.db.execute("UPDATE tasks SET state='pending' WHERE state='blocked'")
                self._refresh(time.time())
            self.db.execute('UPDATE resource_changes SET changed_tasks=? WHERE rowid=?', (len(keys), change))
        return keys

    def snapshot(self):
        with self.lock, self.db:
            self._refresh(time.time())
            return dict(self.db.execute('SELECT state,count(*) FROM tasks GROUP BY state'))

    def results(self, keys=None):
        with self.lock:
            rows = (self.db.execute('SELECT key,state,result FROM tasks ORDER BY rowid') if keys is None
                    else (self.db.execute('SELECT key,state,result FROM tasks WHERE key=?', (key,)).fetchone()
                          for key in keys))
            return {k: dict(state=s, result=_loads_checkpoint(r) if r else None) for k, s, r in rows}

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
        # A separate read connection lets workers renew leases while a large
        # checkpoint is copied. Holding the dispatch lock here can expire
        # healthy workers before their queued heartbeats are accepted.
        with tempfile.TemporaryDirectory(dir=path.parent) as temporary:
            snapshot = Path(temporary)/'queue.sqlite'
            source = sqlite3.connect(f'file:{self.path}?mode=ro', uri=True)
            destination = sqlite3.connect(snapshot)
            try:
                source.backup(destination, pages=256, sleep=.01)
                # The backup inherits WAL mode. Consolidate it before publishing
                # so copying this single file preserves every committed row.
                destination.execute('PRAGMA journal_mode=DELETE')
            finally:
                destination.close()
                source.close()
            os.replace(snapshot, path)

    def close(self):
        with self.lock:
            self.db.close()
