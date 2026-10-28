"""Distribute per-case saved-output metrics, preserving every expected key."""
import argparse
from dataclasses import asdict, is_dataclass
import json
import os
from pathlib import Path
import pickle
import shutil
import sqlite3
import tempfile

from mappa.execution import Task, TaskQueue, read_blob
from .decoder_campaign import save, serve_queue
from .decoder_inventory import signature
from .evidence import sha


def blob_record(blob):
    return asdict(blob) if is_dataclass(blob) else blob


def case_tasks(inventory, manifest, matching, *, memory_bytes, case_seconds):
    matching=Path(matching)
    for case in manifest['inputs']:
        if case['status']!='prepared':
            raise ValueError('Unprepared case needs an explicit metric record, not exclusion')
        index=case['index']
        directions={str(n):r['direction'] for n,r in enumerate(manifest['tasks']) if r['index']==index}
        if set(directions.values())!={'R_to_P','P_to_R'} or len(directions)!=2:
            raise ValueError('Case lacks the two expected directions')
        expected,outputs=[],[]
        for direction in directions:
            for key,state,data in inventory.execute('SELECT key,state,resource FROM tasks WHERE direction=? ORDER BY cut,seed',(int(direction),)):
                expected.append(key)
                if state=='complete':
                    outputs.append(dict(key=key,blob=blob_record(pickle.loads(data)['blob'])))
        if not expected:
            raise ValueError('Case is absent from the expected inventory')
        payload=dict(index=index,expected=expected,outputs=outputs,directions=directions,
            input=dict(path=str(matching/'inputs'/str(index)/'input.json'),sha256=case['input_sha256']),
            reference=dict(path=str(matching/'inputs'/str(index)/'reference.json'),sha256=case['reference_sha256']))
        yield Task(f'metrics/{index:05d}','rebenchmark.decoder_metrics:case_metrics',payload,
                   case=str(index),case_seconds=case_seconds,memory_bytes=memory_bytes)


def prepare(inventory,manifest,matching,queue,reused,memory_bytes,case_seconds,selected=None):
    total=scheduled=0;seen=set();expected_rows=0;accepted_rows=0
    for task in case_tasks(inventory,manifest,matching,memory_bytes=memory_bytes,case_seconds=case_seconds):
        if selected is not None and task.payload['index'] not in selected:
            continue
        total+=1;expected_rows+=len(task.payload['expected']);accepted_rows+=len(task.payload['outputs'])
        if task.key in reused:
            if signature(task)!=signature(reused[task.key]['task']):
                raise ValueError('Reused metric task differs from final input inventory')
            seen.add(task.key)
        else:
            queue.add([task]);scheduled+=1
    if seen!=set(reused):
        raise ValueError('Reused metric result is outside the selected population')
    if selected is not None and total!=len(selected):
        raise ValueError('Selected case is absent')
    if selected is None and expected_rows!=inventory.execute('SELECT count(*) FROM tasks').fetchone()[0]:
        raise ValueError('Case tasks do not cover the expected grid')
    return dict(cases=total,scheduled=scheduled,reused=sorted(seen),expected_rows=expected_rows,
                accepted_output_rows=accepted_rows,missing_output_rows=expected_rows-accepted_rows,
                scientific_audit_accepted=False)


def execute(out):
    plan=json.loads((out/'plan.json').read_text())
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        scratch=Path(scratch);local=scratch/'inventory.sqlite'
        shutil.copyfile(plan['inventory'],local)
        if sha(local)!=plan['inventory_sha256']:
            raise ValueError('Final inventory checksum differs')
        inventory=sqlite3.connect('file:'+str(local)+'?mode=ro&immutable=1',uri=True)
        matching=Path(plan['matching']);manifest_path=matching/'manifest.json'
        if sha(manifest_path)!=plan['manifest_sha256']:
            raise ValueError('Matching manifest checksum differs')
        manifest=json.loads(manifest_path.read_text())
        reused={}
        for path in plan['reuse']:
            source=sqlite3.connect('file:'+str(path)+'?mode=ro&immutable=1',uri=True)
            for key,spec,state,data in source.execute('SELECT key,specification,state,result FROM tasks'):
                if state!='complete' or key in reused:
                    raise ValueError('Incomplete or duplicate prior metric task')
                resource=pickle.loads(data);value=read_blob(resource['blob'],scratch);task=pickle.loads(spec)
                if value['index']!=task.payload['index']:
                    raise ValueError('Prior metric output identity differs')
                reused[key]=dict(task=task,resource=resource)
            source.close()
        queue=TaskQueue(scratch/'queue.sqlite',plan['identity'])
        summary=prepare(inventory,manifest,matching,queue,reused,plan['memory_bytes'],plan['case_seconds'],
                        set(plan['selected']) if 'selected' in plan else None)
        inventory.close();queue.optimize()
        save(out/'selection.json',summary)
        save(out/'reused.json',{key:dict(specification_signature=signature(row['task']),
            blob=blob_record(row['resource']['blob'])) for key,row in reused.items()})
        serve_queue(queue,out,summary);queue.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path)
    execute(parser.parse_args().output)
