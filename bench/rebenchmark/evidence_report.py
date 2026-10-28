"""Join exact saved-output audits to the complete decoder task inventory."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import pickle
import shutil
import sqlite3
import tempfile

from mappa.execution import Blob, Task
from mappa.execution.bundles import read_blobs
from .decoder_inventory import signature
from .decoder_resume import digest
from .evidence_campaign import blob_record


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


class AuditIndex:
    """Missing decoder outputs and unresolved audits stay in the denominator."""
    def __init__(self, path, inventory, selected=None):
        if Path(path).exists():
            raise ValueError('Audit index destination already exists')
        self.db = sqlite3.connect(path)
        self.db.executescript('''
            CREATE TABLE audits (key TEXT PRIMARY KEY, decoder_state TEXT,
                signature TEXT, decoder_blob TEXT, state TEXT, audit_signature TEXT,
                path TEXT, offset INTEGER, size INTEGER, hash TEXT, details TEXT, origin TEXT);
            CREATE INDEX audit_bundle ON audits(path,offset);
            CREATE TABLE sources (origin TEXT PRIMARY KEY, rows INTEGER);
        ''')
        self.population = inventory.execute('SELECT count(*) FROM tasks').fetchone()[0]
        self.selected = selected
        with self.db:
            for key,state,sig,data in inventory.execute('SELECT key,state,signature,resource FROM tasks'):
                if selected is not None and key not in selected:
                    continue
                blob = blob_record(pickle.loads(data)['blob']) if state == 'complete' else None
                self.db.execute('INSERT INTO audits VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                    (key,state,sig,encoded(blob),'missing' if blob else 'no_decoder_output',
                     None,None,None,None,None,None,None))
        count = self.db.execute('SELECT count(*) FROM audits').fetchone()[0]
        if selected is not None and count != len(selected):
            raise ValueError('Selected keys are outside the expected inventory')

    def record(self, task, state, resource, origin, allowed):
        p = task.payload
        key = p['key']
        if task.key != 'audit/'+key or task.operation != 'rebenchmark.decoder_evidence:audit_output':
            raise ValueError('Unexpected audit task identity')
        original = Task(**p['specification'])
        sig = signature(original)
        if (original.key != key or original.operation != 'rebenchmark.decoding:decode_seed'
                or sig != p['specification_signature']):
            raise ValueError('Audit original specification differs')
        previous = self.db.execute(
            'SELECT decoder_state,signature,decoder_blob,state,audit_signature FROM audits WHERE key=?',
            (key,)).fetchone()
        if previous is None or previous[0] != 'complete' or previous[3] not in allowed:
            raise ValueError('Duplicate, unexpected or unavailable decoder output: '+key)
        if previous[1] is not None and previous[1] != sig:
            raise ValueError('Decoder scientific specification differs: '+key)
        if previous[2] != encoded(p['blob']):
            raise ValueError('Audit targets a different decoder output: '+key)
        audit_sig = signature(task)
        if previous[4] is not None and previous[4] != audit_sig:
            raise ValueError('Audit specification changed: '+key)
        blob = Blob(**blob_record(resource['blob'])) if state == 'complete' else None
        self.db.execute('''UPDATE audits SET signature=?,state=?,audit_signature=?,
            path=?,offset=?,size=?,hash=?,origin=?,details=NULL WHERE key=?''',
            (sig,'unchecked' if blob else state,audit_sig,
             blob.path if blob else None,blob.offset if blob else None,
             blob.size if blob else None,blob.sha256 if blob else None,origin,key))

    def queue(self, source, origin, allowed):
        count = 0
        with self.db:
            for key,spec,state,data in source.execute('SELECT key,specification,state,result FROM tasks'):
                task = pickle.loads(spec)
                if key != task.key:
                    raise ValueError('Stored audit task key differs')
                if self.selected is not None and task.payload['key'] not in self.selected:
                    continue
                self.record(task,state,pickle.loads(data) if data else None,origin,allowed)
                count += 1
            self.db.execute('INSERT INTO sources VALUES (?,?)',(origin,count))

    def verify(self, scratch):
        records = self.db.execute('''SELECT path,offset,size,hash FROM audits
            WHERE state='unchecked' ORDER BY path,offset''')
        blobs = (Blob(*row) for row in records)
        count = 0
        for blob,value in read_blobs(blobs,scratch):
            key = value['key']
            row = self.db.execute('''SELECT signature,decoder_blob,state,path,offset,size,hash
                FROM audits WHERE key=?''',(key,)).fetchone()
            if (row is None or row[2] != 'unchecked' or Blob(*row[3:]) != blob
                    or row[0] != value['specification_signature'] or row[1] != encoded(value['blob'])):
                raise ValueError('Saved audit output identity differs: '+key)
            unknown = sum(r['status']=='unknown' for r in value['proof_checks'])
            if (not value['catalogue_verified'] or not value['prefix_bound_verified']
                    or unknown != value['unresolved_proof_checks'] or value['verified'] != (unknown==0)):
                raise ValueError('Saved audit acceptance flags differ: '+key)
            details = {k:value[k] for k in ('verified','candidates','families','decoded_complete',
                'minimum','minimum_proven','incomplete_families','unresolved_proof_checks')}
            self.db.execute('UPDATE audits SET state=?,details=? WHERE key=?',
                ('verified' if value['verified'] else 'unresolved',encoded(details),key))
            count += 1
            if count % 1000 == 0:
                self.db.commit()
        self.db.commit()
        return count

    def summary(self):
        states = dict(self.db.execute('SELECT state,count(*) FROM audits GROUP BY state'))
        expected = self.db.execute('SELECT count(*) FROM audits').fetchone()[0]
        return dict(expected=expected,original_population=self.population,pilot_only=self.selected is not None,
            states=states,sources=dict(self.db.execute('SELECT origin,rows FROM sources')),
            scientific_audit_complete=self.selected is None and states.get('verified',0)==expected,
            scope='Exact accepted-output audit coverage; verified incomplete windows remain incomplete.')


def accept_metrics(value, task, index, smaller_direction):
    """A mode is accepted only when every expected output has a verified audit."""
    p = task.payload
    if (task.operation != 'rebenchmark.decoder_metrics:case_metrics' or value['index'] != p['index']
            or value['input'] != p['input'] or value['reference'] != p['reference']):
        raise ValueError('Metric source identity differs')
    outputs = {r['key']:r['blob'] for r in p['outputs']}
    records = {r['key']:r for r in value['records']}
    if (len(outputs) != len(p['outputs']) or len(set(p['expected'])) != len(p['expected'])
            or len(records) != len(value['records']) or set(records) != set(outputs)):
        raise ValueError('Metric output grid differs')
    if smaller_direction not in ('R_to_P','P_to_R'):
        raise ValueError('Smaller-first direction must come from the verified input sizes')
    accepted = set()
    for key in p['expected']:
        row = index.execute('SELECT decoder_blob,state,details FROM audits WHERE key=?',(key,)).fetchone()
        if row is None:
            raise ValueError('Expected metric key is absent from audit inventory')
        if key in records:
            record = records[key]
            if row[0] != encoded(outputs[key]) or record['blob'] != outputs[key]:
                raise ValueError('Metric and audit decoder outputs differ')
            if row[1] in ('verified','unresolved'):
                details = json.loads(row[2])
                if (record['decoded_complete'] != details['decoded_complete']
                        or record['incomplete_families'] != details['incomplete_families']):
                    raise ValueError('Metric completeness differs from audited evidence')
            if row[1] == 'verified':
                accepted.add(key)
        elif json.loads(row[0]) is not None:
            raise ValueError('Metric task omitted an available decoder output')
    modes = {}
    for name,mode in value['modes'].items():
        scope,direction,prefix = name.split('/')
        if direction == 'smaller_first':
            direction = smaller_direction
        seeds = int(prefix.removeprefix('seeds'))
        def included(key):
            d,_,cut,_,seed = key.split('/')
            return (direction=='bidirectional' or p['directions'][str(int(d))]==direction) and int(seed)<seeds and (scope=='sweep' or int(cut)==0)
        keys = {key for key in p['expected'] if included(key)}
        if mode['expected_rows'] != len(keys):
            raise ValueError('Mode denominator differs')
        modes[name] = dict(mode,scientific_audit_accepted=keys<=accepted,
                          unverified_output_keys=sorted(keys-accepted))
    return dict(index=value['index'],modes=modes,accepted_outputs=len(accepted),
                expected_outputs=len(p['expected']))


def execute(plan, destination):
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        scratch = Path(scratch)
        def snapshot(name, item):
            local = scratch/(name+'.sqlite'); shutil.copyfile(item['path'],local)
            if digest(local) != item['sha256']:
                raise ValueError('Snapshot checksum differs: '+name)
            return sqlite3.connect('file:'+str(local)+'?mode=ro&immutable=1',uri=True)
        source = snapshot('inventory',plan['inventory'])
        path = scratch/'audits.sqlite'
        index = AuditIndex(path,source,set(plan['selected']) if 'selected' in plan else None)
        source.close()
        for n,item in enumerate(plan['queues']):
            source = snapshot('queue-'+str(n),item)
            index.queue(source,item['origin'],set(item['allowed'])); source.close()
        index.verify(scratch)
        summary = index.summary(); index.db.close()
        destination.mkdir(exist_ok=False)
        shutil.copyfile(path,destination/'audits.sqlite')
        summary['sha256'] = digest(destination/'audits.sqlite')
        (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        (destination/'sources.json').write_text(json.dumps(plan,indent=2)+'\n')
        print(json.dumps(summary))


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('plan',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();execute(json.loads(args.plan.read_text()),args.output)
