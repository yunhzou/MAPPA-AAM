"""Retirement must never turn completed or attempted decoder work into reruns."""
from dataclasses import replace
import pickle
import sqlite3
import pytest
from mappa.execution import Task, TaskQueue
from rebenchmark.decoder_resume import prepare_pending, require_retired


def task(key,case='x',seconds=100):
    return Task(key,'rebenchmark.decoding:decode_seed',{'graph':'immutable','seed':3},case=case,case_seconds=seconds)


def test_resume_only_unattempted_and_carry_budget(tmp_path):
    old=TaskQueue(tmp_path/'old.sqlite','old')
    tasks=[task('done'),task('running'),task('lost'),task('error'),task('limited'),task('pending'),task('fresh','new')]
    old.add(tasks)
    for t,state in zip(tasks,['complete','running','lost','error','timeout']):
        old.db.execute('UPDATE tasks SET state=? WHERE key=?',(state,t.key))
        old.db.execute('INSERT INTO attempts VALUES (?,?,?,?,?,?,?)',(t.key,t.key,'owner',10,20,state,None))
    old.db.execute("UPDATE cases SET deadline=100 WHERE name='x'");old.db.commit()
    new=TaskQueue(tmp_path/'new.sqlite','new')
    summary=prepare_pending(old.db,new,40)
    assert summary['scheduled']==2
    rows={k:pickle.loads(s) for k,s in new.db.execute('SELECT key,specification FROM tasks')}
    assert rows=={'pending':replace(tasks[5],case_seconds=60),'fresh':tasks[6]}
    assert dict(old.db.execute('SELECT key,state FROM tasks'))['done']=='complete'
    assert not list(new.db.execute('SELECT * FROM attempts'))
    assert new.db.execute("SELECT deadline FROM cases WHERE name='x'").fetchone()==(None,)


def test_exhausted_direction_is_not_reset(tmp_path):
    old=TaskQueue(tmp_path/'old.sqlite','old');old.add([task('a'),task('b','fresh')])
    old.db.execute("UPDATE cases SET deadline=50 WHERE name='x'");old.db.commit()
    new=TaskQueue(tmp_path/'new.sqlite','new');r=prepare_pending(old.db,new,60)
    assert r['expired_pending']=={'x':1} and r['scheduled']==1
    assert list(new.db.execute('SELECT key FROM tasks'))==[('b',)]


def test_pending_with_prior_attempt_rejected(tmp_path):
    old=TaskQueue(tmp_path/'old.sqlite','old');old.add([task('a')])
    old.db.execute('INSERT INTO attempts VALUES (?,?,?,?,?,?,?)',('token','a','owner',10,20,'lost',None));old.db.commit()
    new=TaskQueue(tmp_path/'new.sqlite','new')
    with pytest.raises(ValueError,match='prior ownership'):prepare_pending(old.db,new,30)


@pytest.mark.parametrize('coordinator,worker,live,allowed',[
    ('RUNNING','COMPLETED','RUNNING',False),('TIMEOUT','RUNNING','RUNNING',False),
    ('TIMEOUT','COMPLETING','COMPLETING',False),('TIMEOUT','PENDING','PENDING',True),
    ('COMPLETED','COMPLETED','',True),('TIMEOUT','FAILED','RUNNING',False)])
def test_original_owners_must_retire(monkeypatch,coordinator,worker,live,allowed):
    from types import SimpleNamespace
    def run(cmd,**kwargs):
        return SimpleNamespace(stdout=f'1|{coordinator}|end\n2_0|{worker}|end\n' if cmd[0]=='sacct' else live)
    monkeypatch.setattr('rebenchmark.decoder_resume.subprocess.run',run)
    if allowed:assert require_retired(1,2)['checked_at']>0
    else:
        with pytest.raises(ValueError):require_retired(1,2)
