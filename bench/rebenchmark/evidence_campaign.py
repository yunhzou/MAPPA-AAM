"""Schedule saved-output verification from the final, complete task inventory."""
import argparse
from dataclasses import asdict, is_dataclass
import json
import os
from pathlib import Path
import pickle
import shutil
import sqlite3
import tempfile

from mappa.execution import Blob, Task, TaskQueue, read_blob
from .decoder_campaign import save, serve_queue
from .decoder_inventory import signature
from .decoder_resume import digest


def blob_record(blob):
    return asdict(blob) if is_dataclass(blob) else blob


def prepare(inventory, original, pilots, reused, queue, proof_seconds, memory_bytes,
            case_seconds, selected=None):
    """Keep unresolved executions in the denominator and never repeat audits."""
    expected = scheduled = 0
    missing, accepted_reuse, batch = [], {}, []
    for key, state, sig, data in inventory.execute(
            'SELECT key,state,signature,resource FROM tasks ORDER BY direction,cut,seed'):
        if selected is not None and key not in selected:
            continue
        expected += 1
        if state != 'complete':
            missing.append(dict(key=key, state=state)); continue
        record = original.execute('SELECT specification FROM tasks WHERE key=?', (key,)).fetchone()
        task = pickle.loads(record[0]) if record else pilots[key]
        current = signature(task)
        if task.key != key or task.operation != 'rebenchmark.decoding:decode_seed' or (sig and current != sig):
            raise ValueError('Original task specification differs: '+key)
        blob = blob_record(pickle.loads(data)['blob'])
        if key in reused:
            row = reused[key]
            if (not row['verified'] or row['blob'] != blob or
                    row['specification_signature'] != current):
                raise ValueError('Reused audit differs from accepted output: '+key)
            accepted_reuse[key] = row
            continue
        payload = dict(key=key, blob=blob, specification=asdict(task),
                       specification_signature=current, proof_seconds=proof_seconds)
        batch.append(Task('audit/'+key, 'rebenchmark.decoder_evidence:audit_output', payload,
                          case=key.split('/')[0], case_seconds=case_seconds, memory_bytes=memory_bytes))
        scheduled += 1
        if len(batch) == 1000:
            queue.add(batch); batch = []
    if batch:
        queue.add(batch)
    if selected is not None and expected != len(selected):
        raise ValueError('Selected audit keys are absent from the inventory')
    if set(accepted_reuse) != set(reused):
        raise ValueError('Reused audit is absent from this inventory selection')
    return dict(expected=expected, scheduled=scheduled, missing=missing,
                reused=accepted_reuse, scientific_audit_complete=False)


def execute(out):
    plan = json.loads((out/'plan.json').read_text())
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        scratch = Path(scratch)
        sources = []
        def snapshot(name, source, expected=None):
            dest = scratch/(name+'.sqlite')
            if expected:
                shutil.copyfile(source, dest)
                if digest(dest) != expected:
                    raise ValueError('Source snapshot checksum differs: '+name)
            else:
                # Old pilot queues may retain WAL journals. Copy committed state.
                db = sqlite3.connect('file:'+str(source)+'?mode=ro', uri=True)
                target = sqlite3.connect(dest); db.backup(target); target.close(); db.close()
            db = sqlite3.connect('file:'+str(dest)+'?mode=ro&immutable=1', uri=True)
            if db.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
                raise ValueError('Invalid source database: '+name)
            sources.append(dict(name=name, source=str(source), sha256=digest(dest)))
            return db
        inventory = snapshot('inventory', plan['inventory'], plan['inventory_sha256'])
        original = snapshot('original', plan['original'], plan['original_sha256'])
        pilots = {}
        for number, path in enumerate(plan['pilots']):
            db = snapshot('pilot-'+str(number), path)
            for key, spec in db.execute('SELECT key,specification FROM tasks'):
                if key in pilots:
                    raise ValueError('Duplicate pilot task')
                pilots[key] = pickle.loads(spec)
            db.close()
        reused = {}
        for path in plan['audit_reuse']:
            for row in json.loads(Path(path).read_text())['rows']:
                key = row['key']
                if key in reused:
                    raise ValueError('Duplicate prior audit')
                saved = read_blob(Blob(**row['audit_blob']), scratch)
                if any(saved[k] != row[k] for k in ('key','blob','specification_signature','verified')):
                    raise ValueError('Prior audit receipt differs from verified bundle')
                reused[key] = row
        queue = TaskQueue(scratch/'queue.sqlite', plan['identity'])
        summary = prepare(inventory, original, pilots, reused, queue,
                          plan['proof_seconds'], plan['memory_bytes'], plan['case_seconds'],
                          set(plan['selected']) if 'selected' in plan else None)
        inventory.close(); original.close()
        save(out/'sources.json', sources)
        save(out/'selection.json', summary)
        queue.optimize()
        serve_queue(queue, out, {k:summary[k] for k in ('expected','scheduled')})
        queue.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('output', type=Path)
    execute(parser.parse_args().output)
