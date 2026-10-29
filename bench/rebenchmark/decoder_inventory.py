"""Reconcile decoder execution records before auditing scientific outputs.

The expected cut/seed grid owns the denominator. Continuations and repairs
replace only their diagnosed predecessor states; completed outputs cannot be
silently replaced. This inventory does not certify event-window completeness.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import pickle
import sqlite3
from mappa.artifacts import _loads_checkpoint


def signature(task):
    record = asdict(task)
    record.pop('case_seconds')  # Continuations retain a smaller remaining budget.
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()


class Inventory:
    def __init__(self, path, directions, seeds=10):
        if Path(path).exists():
            raise ValueError('Inventory destination already exists')
        self.db = sqlite3.connect(path)
        self.db.executescript('''
            CREATE TABLE tasks (key TEXT PRIMARY KEY, direction INTEGER, cut INTEGER,
                seed INTEGER, state TEXT NOT NULL, signature TEXT, resource BLOB, origin TEXT);
            CREATE INDEX direction_tasks ON tasks(direction);
            CREATE TABLE sources (origin TEXT PRIMARY KEY, rows INTEGER);
        ''')
        with self.db:
            for row in directions:
                if row['expected'] != row['cuts'] * seeds:
                    raise ValueError('Direction does not contain the full cut/seed grid')
                self.db.executemany('INSERT INTO tasks VALUES (?,?,?,?,?,?,?,?)',
                    ((f'{row["task"]:05d}/cut/{cut:05d}/seed/{seed:02d}',
                      row['task'], cut, seed, 'missing', None, None, None)
                     for cut in range(row['cuts']) for seed in range(seeds)))

    def replace(self, key, state, resource, origin, allowed, task=None):
        previous = self.db.execute('SELECT state,signature FROM tasks WHERE key=?', (key,)).fetchone()
        if previous is None or previous[0] not in allowed:
            raise ValueError('Unexpected replacement or duplicate output: '+key)
        digest = previous[1]
        if task is not None:
            p = task.payload
            expected = f'{p["task_id"]:05d}/cut/{p["cut"]:05d}/seed/{p["seed"]:02d}'
            if task.key != key or key != expected or task.operation != 'rebenchmark.decoding:decode_seed':
                raise ValueError('Task identity differs: '+key)
            current = signature(task)
            if digest is not None and current != digest:
                raise ValueError('Scientific task specification changed: '+key)
            digest = current
        if state == 'complete' and (not resource or 'blob' not in resource):
            raise ValueError('Completed execution lacks its output: '+key)
        self.db.execute('UPDATE tasks SET state=?,signature=?,resource=?,origin=? WHERE key=?',
            (state, digest, pickle.dumps(resource, protocol=5), origin, key))

    def queue(self, path, origin, allowed):
        source = sqlite3.connect('file:'+str(Path(path).resolve())+'?mode=ro&immutable=1', uri=True)
        count = 0
        with self.db:
            for key, spec, state, result in source.execute('SELECT key,specification,state,result FROM tasks'):
                self.replace(key, state, _loads_checkpoint(result) if result else None,
                             origin, allowed, _loads_checkpoint(spec))
                count += 1
            self.db.execute('INSERT INTO sources VALUES (?,?)', (origin, count))
        source.close()

    def recovered(self, rows, origin):
        with self.db:
            for key, row in rows.items():
                if row['key'] != key:
                    raise ValueError('Recovery key differs')
                self.replace(key, 'complete', {'blob': row['blob']}, origin, {'lost'})
            self.db.execute('INSERT INTO sources VALUES (?,?)', (origin, len(rows)))

    def pilots(self, rows, origin):
        with self.db:
            for key, row in rows.items():
                if row['state'] != 'complete':
                    raise ValueError('Pilot execution is not complete')
                self.replace(key, 'complete', row['resource'], origin, {'missing'})
            self.db.execute('INSERT INTO sources VALUES (?,?)', (origin, len(rows)))

    def summary(self):
        return dict(expected=self.db.execute('SELECT count(*) FROM tasks').fetchone()[0],
            states=dict(self.db.execute('SELECT state,count(*) FROM tasks GROUP BY state')),
            sources=dict(self.db.execute('SELECT origin,rows FROM sources')),
            directions=[dict(direction=d, state=s, tasks=n) for d,s,n in self.db.execute(
                'SELECT direction,state,count(*) FROM tasks GROUP BY direction,state ORDER BY direction,state')],
            scientific_audit_complete=False,
            scope='Execution accounting only. Every published output, witness and completeness proof still needs scientific audit.')


def execute(plan, output):
    """All queue paths must be standalone snapshots, never live WAL files."""
    definition = json.loads(Path(plan['expected']).read_text())
    inventory = Inventory(output, definition['directions'], plan.get('seeds', 10))
    inventory.queue(plan['original'], 'original', {'missing'})
    inventory.pilots(json.loads(Path(plan['pilots']).read_text()), 'pilots')
    recovery = json.loads(Path(plan['recovery']).read_text())
    if recovery['missing'] or recovery['duplicates']:
        raise ValueError('Recovery receipt has unresolved collisions or missing targets')
    inventory.recovered(recovery['recovered'], 'startup_recovery')
    inventory.queue(plan['continuation'], 'continuation', {'pending'})
    for repair in plan['repairs']:
        selection = json.loads(Path(repair['selection']).read_text())
        expected = {r['key'] for r in selection['selected']}
        source = sqlite3.connect('file:'+str(Path(repair['queue']).resolve())+'?mode=ro&immutable=1', uri=True)
        actual = {r[0] for r in source.execute('SELECT key FROM tasks')}
        source.close()
        if actual != expected:
            raise ValueError('Repair queue differs from exact selection')
        inventory.queue(repair['queue'], repair['name'], {'error'})
    summary = inventory.summary()
    if summary['expected'] != definition['total'] or summary['states'].get('missing', 0):
        raise ValueError('Expected task coverage is incomplete')
    inventory.db.close()
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('plan', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    summary = execute(json.loads(args.plan.read_text()), args.output)
    args.output.with_suffix('.json').write_text(json.dumps(summary, indent=2)+'\n')
