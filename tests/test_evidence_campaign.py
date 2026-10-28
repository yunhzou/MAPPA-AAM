from dataclasses import asdict, replace
import pickle
import pytest
from mappa.execution import TaskQueue
from rebenchmark.decoder_inventory import Inventory, signature
from rebenchmark.evidence_campaign import prepare
from test_decoder_inventory import task


def setup(tmp_path):
    inventory = Inventory(tmp_path/'i.sqlite', [dict(task=0,cuts=1,expected=4)], seeds=4)
    original = TaskQueue(tmp_path/'o.sqlite', 'original')
    original.add([task(0), task(1), task(2)])
    for seed,state in enumerate(['complete','error','complete','complete']):
        inventory.replace(task(seed).key,state,dict(blob=dict(path='saved',offset=seed,size=1,sha256='hash')),
                          'source',{'missing'},task(seed) if seed!=3 else None)
    queue = TaskQueue(tmp_path/'q.sqlite', 'audit')
    return inventory,original,queue


def test_reuse_partial_denominator_pilot_and_original_metadata(tmp_path):
    i,o,q=setup(tmp_path)
    prior=dict(key=task(0).key,verified=True,blob=dict(path='saved',offset=0,size=1,sha256='hash'),
               specification_signature=signature(task(0)))
    result=prepare(i.db,o.db,{task(3).key:task(3)},{task(0).key:prior},q,10,4096,3600)
    assert result['expected']==4 and result['scheduled']==2
    assert result['missing']==[dict(key=task(1).key,state='error')]
    specs=[pickle.loads(r[0]) for r in q.db.execute('SELECT specification FROM tasks')]
    assert [s.payload['specification'] for s in specs]==[asdict(task(2)),asdict(task(3))]
    assert not result['scientific_audit_complete']


@pytest.mark.parametrize('change',['blob','signature','unverified','missing'])
def test_reject_reuse_that_is_not_the_same_verified_output(tmp_path,change):
    i,o,q=setup(tmp_path)
    prior=dict(key=task(0).key,verified=True,blob=dict(path='saved',offset=0,size=1,sha256='hash'),
               specification_signature=signature(task(0)))
    if change=='blob':prior['blob']['offset']=99
    if change=='signature':prior['specification_signature']='wrong'
    if change=='unverified':prior['verified']=False
    reused={task(9).key if change=='missing' else task(0).key:prior}
    with pytest.raises(ValueError,match='audit'):
        prepare(i.db,o.db,{task(3).key:task(3)},reused,q,10,4096,3600)


def test_changed_task_and_absent_selection_rejected(tmp_path):
    i,o,q=setup(tmp_path)
    with pytest.raises(ValueError,match='Selected'):
        prepare(i.db,o.db,{}, {},q,10,4096,3600,{task(9).key})
    with o.db:
        o.db.execute('UPDATE tasks SET specification=? WHERE key=?',
            (pickle.dumps(replace(task(0),payload=dict(task(0).payload,slack=6))),task(0).key))
    with pytest.raises(ValueError,match='specification'):
        prepare(i.db,o.db,{}, {},q,10,4096,3600,{task(0).key})
