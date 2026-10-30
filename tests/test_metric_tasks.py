from dataclasses import replace
import pytest
from mappa.mappa_core.core.scheduler import TaskQueue
from rebenchmark.decoder_inventory import Inventory
from rebenchmark.metric_tasks import case_tasks,prepare
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


def test_fresh_grid_preserves_timeout_and_rejects_changed_task(tmp_path):
    import pickle
    from dataclasses import asdict
    from mappa.mappa_core.core.scheduler import Task, BundleWriter
    from rebenchmark.metric_tasks import fresh_case_tasks
    writer=BundleWriter(tmp_path/'bundles',tmp_path)
    source=TaskQueue(tmp_path/'fresh.sqlite','fresh')
    directions=[]
    for direction in (0,1):
        key=f'{direction:05d}/cut/00000'
        parent=Task(key,'mappa.mappa_core.core.scheduler.search:finalize_cut',{},case=f'{direction:05d}')
        child=Task(key+'/seed/00','rebenchmark.operations:convert_seed',{'cut':key},case=f'{direction:05d}')
        source.add([parent,child])
        resource=pickle.dumps({'blob':writer.put({'saved':True})})
        source.db.execute("UPDATE tasks SET state='complete',result=? WHERE case_name=?",(resource,f'{direction:05d}'))
        if direction:
            source.db.execute("UPDATE tasks SET state='timeout',result=NULL WHERE key=?",(child.key,))
        directions.append(dict(direction=direction,case=0,reversed=bool(direction),blob=asdict(writer.put([parent,child]))))
    source.db.commit()
    manifest={'inputs':[dict(index=0,input_sha256='input',reference_sha256='ref')]}
    index=dict(tasks=4,directions=directions)
    tasks=list(fresh_case_tasks(source.db,index,manifest,tmp_path,tmp_path,memory_bytes=4096,case_seconds=3600))
    assert len(tasks)==1
    assert len(tasks[0].payload['expected'])==2
    assert len(tasks[0].payload['outputs'])==1
    assert tasks[0].payload['missing_tasks'][0]['state']=='timeout'
    changed=replace(child,payload={'cut':'wrong'})
    source.db.execute('UPDATE tasks SET specification=? WHERE key=?',(pickle.dumps(changed),child.key))
    with pytest.raises(ValueError,match='scientific task changed'):
        list(fresh_case_tasks(source.db,index,manifest,tmp_path,tmp_path,memory_bytes=4096,case_seconds=3600))
