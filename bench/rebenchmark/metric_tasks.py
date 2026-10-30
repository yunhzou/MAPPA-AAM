"""Distribute per-case saved-output metrics, preserving every expected key."""
from dataclasses import asdict, is_dataclass
from pathlib import Path

from mappa.mappa_core.core.scheduler import Task, read_blob
from .decoder_inventory import signature
from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint


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
                    outputs.append(dict(key=key,blob=blob_record(_loads_checkpoint(data)['blob'])))
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


def fresh_case_tasks(source, index, manifest, matching, scratch, *, memory_bytes, case_seconds, selected=None):
    """Reconcile the frozen fresh task grid and retain every missing output."""
    directions = {row['direction']: row for row in index['directions']}
    if source.execute('SELECT count(*) FROM tasks').fetchone()[0] != index['tasks']:
        raise ValueError('Fresh queue has a different task denominator')
    from mappa.mappa_core.core.scheduler import Blob
    for case in manifest['inputs']:
        number = case['index']
        if selected is not None and number not in selected:
            continue
        expected, outputs, search, missing = [], [], [], []
        for direction in (2*number, 2*number+1):
            row = directions[direction]
            if row['case'] != number or row['reversed'] != bool(direction % 2):
                raise ValueError('Fresh direction identity differs')
            planned = read_blob(Blob(**row['blob']), scratch)
            saved = {k:(s, state, result) for k,s,state,result in source.execute(
                'SELECT key,specification,state,result FROM tasks WHERE case_name=?', (f'{direction:05d}',))}
            if set(saved) != {t.key for t in planned}:
                raise ValueError('Fresh task keys differ from the immutable index')
            for task in planned:
                specification, state, result = saved[task.key]
                if _loads_checkpoint(specification) != task:
                    raise ValueError('Fresh scientific task changed: '+task.key)
                if state in ('pending', 'running'):
                    raise ValueError('Fresh computation still owns live work')
                conversion = task.operation == 'rebenchmark.operations:convert_seed'
                if conversion:
                    expected.append(task.key)
                if state != 'complete':
                    missing.append(dict(key=task.key,state=state,operation=task.operation))
                    continue
                resource = _loads_checkpoint(result)
                item = dict(key=task.key,blob=blob_record(resource['blob']),
                    resource={k:v for k,v in resource.items() if k!='blob'},payload=task.payload)
                if conversion:
                    parent=saved[task.payload['cut']]
                    if parent[1]!='complete':
                        raise ValueError('Conversion has no accepted finalized parent')
                    item['finalized_blob']=blob_record(_loads_checkpoint(parent[2])['blob'])
                (outputs if conversion else search).append(item)
        folder = Path(matching)/'inputs'/str(number)
        payload = dict(index=number,expected=expected,outputs=outputs,search_outputs=search,
            missing_tasks=missing,source_kind='fresh_conversion',
            directions={str(2*number):'R_to_P',str(2*number+1):'P_to_R'},
            input=dict(path=str(folder/'input.json'),sha256=case['input_sha256']),
            reference=dict(path=str(folder/'reference.json'),sha256=case['reference_sha256']))
        yield Task(f'metrics/{number:05d}','rebenchmark.decoder_metrics:case_metrics',payload,
            case=str(number),case_seconds=case_seconds,memory_bytes=memory_bytes)
