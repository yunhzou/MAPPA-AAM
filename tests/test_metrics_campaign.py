from dataclasses import replace
import pytest
from mappa.mappa_core.core.scheduler import TaskQueue
from rebenchmark.decoder_inventory import Inventory
from rebenchmark.metrics_campaign import case_tasks,prepare
from test_decoder_inventory import task


def setup(tmp_path):
    i=Inventory(tmp_path/'i.sqlite',[dict(task=0,cuts=1,expected=1),dict(task=1,cuts=1,expected=1)],seeds=1)
    t=task();i.replace(t.key,'complete',{'blob':dict(path='saved',offset=0,size=1,sha256='hash')},'o',{'missing'},t)
    missing=replace(t,key='00001/cut/00000/seed/00',payload=dict(t.payload,task_id=1))
    i.replace(missing.key,'error',None,'o',{'missing'},missing)
    m=dict(inputs=[dict(index=0,status='prepared',input_sha256='input',reference_sha256='reference')],
           tasks=[dict(index=0,direction=d) for d in ['R_to_P','P_to_R']])
    q=TaskQueue(tmp_path/'q.sqlite','metrics')
    return i,m,q


def test_case_grid_preserves_missing_direction_outputs_and_reuses(tmp_path):
    i,m,q=setup(tmp_path)
    t=list(case_tasks(i.db,m,'source',memory_bytes=4096,case_seconds=3600))[0]
    assert len(t.payload['expected'])==2 and len(t.payload['outputs'])==1
    summary=prepare(i.db,m,'source',q,{t.key:dict(task=t)},4096,3600)
    assert summary['scheduled']==0 and summary['missing_output_rows']==1
    assert not summary['scientific_audit_accepted']
    with pytest.raises(ValueError,match='differs'):
        prepare(i.db,m,'changed',q,{t.key:dict(task=t)},4096,3600)


def test_unprepared_or_missing_direction_is_not_dropped(tmp_path):
    i,m,q=setup(tmp_path)
    m['inputs'][0]['status']='preparation_error'
    with pytest.raises(ValueError,match='Unprepared'):list(case_tasks(i.db,m,'source',memory_bytes=4096,case_seconds=3600))
    m['inputs'][0]['status']='prepared';m['tasks'].pop()
    with pytest.raises(ValueError,match='directions'):list(case_tasks(i.db,m,'source',memory_bytes=4096,case_seconds=3600))
