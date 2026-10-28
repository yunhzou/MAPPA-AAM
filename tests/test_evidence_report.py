from copy import deepcopy
from dataclasses import asdict,replace
import pickle
import sqlite3

import pytest
from mappa.execution import Blob,BundleWriter,Task
from mappa.execution.bundles import read_blobs
from rebenchmark.decoder_inventory import signature
from rebenchmark.evidence_report import AuditIndex,accept_metrics


def setup(tmp_path, *, missing=False):
    source=sqlite3.connect(':memory:')
    source.execute('CREATE TABLE tasks (key TEXT,state TEXT,signature TEXT,resource BLOB)')
    writer=BundleWriter(tmp_path/'bundles',tmp_path)
    tasks=[]
    for n in range(2):
        original=Task(f'{n:05d}/cut/00000/seed/00','rebenchmark.decoding:decode_seed',dict(task_id=n))
        blob=asdict(writer.put(dict(decoded=n)))
        source.execute('INSERT INTO tasks VALUES (?,?,?,?)',(original.key,'error' if missing and n else 'complete',
            signature(original),pickle.dumps(dict(blob=blob))))
        tasks.append(Task('audit/'+original.key,'rebenchmark.decoder_evidence:audit_output',
            dict(key=original.key,blob=blob,specification=asdict(original),specification_signature=signature(original),proof_seconds=10),case=str(n)))
    return AuditIndex(tmp_path/'index.sqlite',source),tasks,writer


def output(task, *, complete=True, unknown=False):
    return dict(key=task.payload['key'],blob=task.payload['blob'],specification_signature=task.payload['specification_signature'],
        verified=not unknown,candidates=1,families=1,decoded_complete=complete,minimum=2,minimum_proven=complete,
        incomplete_families=[] if complete else ['limited'],proof_checks=[dict(status='unknown' if unknown else 'sat')],
        unresolved_proof_checks=int(unknown),catalogue_verified=True,prefix_bound_verified=True)


def put(index,task,writer,**kwargs):
    blob=writer.put(output(task,**kwargs))
    index.record(task,'complete',dict(blob=blob),'main',{'missing'})
    return blob


def test_full_grid_and_verified_incomplete_are_distinct(tmp_path):
    index,tasks,writer=setup(tmp_path,missing=True)
    put(index,tasks[0],writer,complete=False)
    index.verify(tmp_path)
    s=index.summary()
    assert s['expected']==2 and s['states']=={'no_decoder_output':1,'verified':1}
    assert not s['scientific_audit_complete']
    assert 'false' in index.db.execute('SELECT details FROM audits WHERE state=?',('verified',)).fetchone()[0]


def test_unknown_proof_never_becomes_verified(tmp_path):
    index,tasks,writer=setup(tmp_path)
    put(index,tasks[0],writer,unknown=True);put(index,tasks[1],writer)
    assert index.verify(tmp_path)==2
    assert index.summary()['states']=={'unresolved':1,'verified':1}
    assert not index.summary()['scientific_audit_complete']


def test_repairs_exclude_successes_and_preserve_specifications(tmp_path):
    index,tasks,writer=setup(tmp_path)
    index.record(tasks[0],'error',None,'old',{'missing'})
    fixed=replace(tasks[0],case_seconds=100)
    index.record(fixed,'complete',dict(blob=writer.put(output(fixed))),'repair',{'error'})
    with pytest.raises(ValueError,match='Duplicate'):
        index.record(fixed,'complete',dict(blob=writer.put(output(fixed))),'duplicate',{'error'})
    altered=replace(tasks[1],payload=dict(tasks[1].payload,blob=dict(tasks[1].payload['blob'],sha256='wrong')))
    with pytest.raises(ValueError,match='different decoder output'):
        index.record(altered,'complete',dict(blob=writer.put(output(altered))),'bad',{'missing'})


def test_saved_audit_identity_and_flags_are_checked(tmp_path):
    index,tasks,writer=setup(tmp_path)
    value=output(tasks[0]);value['blob']=tasks[1].payload['blob']
    index.record(tasks[0],'complete',dict(blob=writer.put(value)),'main',{'missing'})
    with pytest.raises(ValueError,match='identity differs'):index.verify(tmp_path)


def test_adjacent_verified_ranges_reuse_one_file_and_reject_corruption(tmp_path,monkeypatch):
    writer=BundleWriter(tmp_path/'bundles',tmp_path)
    blobs=[writer.put('large'*1000),writer.put('small')]
    import builtins
    real=builtins.open;opened=[]
    def counted(path,*args,**kwargs):
        if str(path)==blobs[0].path:opened.append(str(path))
        return real(path,*args,**kwargs)
    monkeypatch.setattr(builtins,'open',counted)
    assert [v for _,v in read_blobs(blobs,tmp_path)]==['large'*1000,'small']
    assert len(opened)==1
    with pytest.raises(ValueError,match='checksum'):
        list(read_blobs([replace(blobs[1],sha256='bad')],tmp_path))


def test_metric_acceptance_uses_exact_blobs_and_each_direction(tmp_path):
    index,tasks,writer=setup(tmp_path,missing=True)
    put(index,tasks[0],writer);index.verify(tmp_path)
    keys=[t.payload['key'] for t in tasks]
    task=Task('metrics/00000','rebenchmark.decoder_metrics:case_metrics',dict(index=0,input={'sha256':'i'},reference={'sha256':'r'},
        expected=keys,outputs=[dict(key=keys[0],blob=tasks[0].payload['blob'])],directions={'0':'R_to_P','1':'P_to_R'}))
    value=dict(index=0,input=task.payload['input'],reference=task.payload['reference'],records=[dict(key=keys[0],blob=tasks[0].payload['blob'],decoded_complete=True,incomplete_families=[])],
        modes={'uncut/R_to_P/seeds1':dict(expected_rows=1),'uncut/P_to_R/seeds1':dict(expected_rows=1),
               'uncut/bidirectional/seeds1':dict(expected_rows=2),'uncut/smaller_first/seeds1':dict(expected_rows=1)})
    result=accept_metrics(value,task,index.db,'P_to_R')
    assert result['modes']['uncut/R_to_P/seeds1']['scientific_audit_accepted']
    assert not result['modes']['uncut/P_to_R/seeds1']['scientific_audit_accepted']
    assert not result['modes']['uncut/bidirectional/seeds1']['scientific_audit_accepted']
    assert not result['modes']['uncut/smaller_first/seeds1']['scientific_audit_accepted']
    changed=deepcopy(value);changed['records'][0]['blob']['sha256']='wrong'
    with pytest.raises(ValueError,match='outputs differ'):accept_metrics(changed,task,index.db,'P_to_R')
