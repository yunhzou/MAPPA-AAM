"""Read timing records from archives without extracting shared-storage files."""
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile


def timing_record(progress):
    state=json.loads(progress.read_text())
    if not state['complete']:
        return dict(task=state['task'],complete=False,stages=state['stages'])
    selected={};hashes={};manifest=None
    with tarfile.open(state['archive']['path'],'r:') as archive:
        for member in archive:
            name=member.name
            if name=='archive-manifest.json':
                manifest=json.load(archive.extractfile(member))['files']
            elif name in ('search.json','search-resources.json','analysis-resources.json',
                          'task.json','environment.json') or (
                    name.startswith('analysis/') and name.endswith('/complete.json')):
                data=archive.extractfile(member).read()
                hashes[name]=hashlib.sha256(data).hexdigest()
                selected[name]=json.loads(data)
    assert manifest is not None
    for name,digest in hashes.items():
        assert manifest[name]['sha256']==digest,(progress,name)
    search=selected['search.json']
    finalization=[v['finalization'] for n,v in selected.items()
                  if n.startswith('analysis/') and n.endswith('/complete.json')]
    assert len(finalization)==search['metrics']['cut_count']
    return dict(task=state['task'],complete=True,index=state['index'],direction=state['direction'],
        search=search['timing'],metrics=search['metrics'],
        search_resources=selected['search-resources.json'],
        analysis_resources=selected['analysis-resources.json'],
        finalization_cpu_seconds=sum(v['cpu_seconds'] for v in finalization),
        finalization_summed_worker_wall_seconds=sum(v['wall_seconds'] for v in finalization),
        environment=selected['environment.json'],archive_sha256=state['archive']['sha256'])


def main(run,output):
    paths=sorted((Path(run)/'progress').glob('*.json'))
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        records=[]
        for record in pool.map(timing_record,paths,chunksize=8):
            records.append(record)
            if len(records)%200==0:
                print(f'{len(records)}/{len(paths)} timing records',flush=True)
    result=dict(records=records,probe_cpu=subprocess.check_output(['lscpu'],text=True),
        source_manifest_sha256=hashlib.sha256((Path(run)/'manifest.json').read_bytes()).hexdigest(),
        scope='Saved timings only; no search or decoding rerun. Four archive readers; no extraction.')
    with gzip.open(output,'wt') as stream:
        json.dump(result,stream,separators=(',',':'))


def distribution(values):
    """Measured samples only; a sum of worker wall samples is not latency."""
    import numpy as np
    values=np.asarray(values,dtype=float)
    if not len(values):
        return dict(count=0,total=0.,median=None,p95=None,maximum=None)
    return dict(count=len(values),total=float(values.sum()),median=float(np.median(values)),
                p95=float(np.percentile(values,95)),maximum=float(values.max()))


def case_timing(record):
    """Original science timings, separate from later verification work."""
    from collections import defaultdict
    phases=defaultdict(lambda:dict(cpu=[],wall=[],peak=[]))
    directions=defaultdict(lambda:defaultdict(float));prefixes=defaultdict(lambda:defaultdict(float))
    def add(name,cpu,wall,peak,direction,cut,seed=None):
        phases[name]['cpu'].append(cpu);phases[name]['wall'].append(wall);phases[name]['peak'].append(peak)
        directions[direction][name+'_cpu']+=cpu
        directions[direction][name+'_worker_wall']+=wall
        for n in range(1,11):
            if seed is not None and seed>=n:continue
            for scope in (('uncut','sweep') if cut==0 else ('sweep',)):
                prefixes[f'{scope}/{direction}/seeds{n}'][name+'_cpu']+=cpu
                prefixes[f'{scope}/{direction}/seeds{n}'][name+'_worker_wall']+=wall
    for row in record['conversion_checks']:
        d,_,cut,_,seed=row['key'].split('/');direction='P_to_R' if int(d)%2 else 'R_to_P'
        resource=row['resource'];peak=resource['peak_tree_rss_bytes']
        for name,timing in row['timing'].items():
            add('conversion_'+name,timing['cpu_seconds'],timing['wall_seconds'],peak,direction,int(cut),int(seed))
        add('conversion_task',resource['cpu_seconds'],resource['compute_seconds'],peak,direction,int(cut),int(seed))
        add('conversion_publication',resource['publish_cpu_seconds'],resource['publish_seconds'],peak,direction,int(cut),int(seed))
        population='complete_conversion' if row['complete'] else 'incomplete_conversion'
        phases[population]['cpu'].append(row['timing']['decoding']['cpu_seconds'])
        phases[population]['wall'].append(row['timing']['decoding']['wall_seconds'])
        phases[population]['peak'].append(peak)
    for row in record['search_resources']:
        parts=row['key'].split('/');direction='P_to_R' if int(parts[0])%2 else 'R_to_P'
        phase='matching_batch' if parts[3]=='seeds' else 'finalization'
        seed=int(parts[4]) if phase=='matching_batch' else None
        add(phase,row['cpu_seconds'],row['compute_seconds'],row['peak_tree_rss_bytes'],direction,int(parts[2]),seed)
        add(phase+'_publication',row['publish_cpu_seconds'],row['publish_seconds'],row['peak_tree_rss_bytes'],direction,int(parts[2]),seed)
    return dict(phases=dict(phases),directions=dict(directions),prefixes=dict(prefixes))


def execution_timeline(queue_path):
    """Claim-to-accept intervals include dispatch and publication, not pure compute."""
    import sqlite3
    import numpy as np
    db=sqlite3.connect('file:'+str(queue_path)+'?mode=ro&immutable=1',uri=True)
    intervals=np.array(db.execute('SELECT started,finished FROM attempts WHERE finished IS NOT NULL').fetchall())
    states=dict(db.execute('SELECT state,count(*) FROM attempts GROUP BY state'));db.close()
    if not len(intervals):return dict(attempts=states)
    points=np.concatenate((intervals[:,0],intervals[:,1]))
    changes=np.concatenate((np.ones(len(intervals),dtype=int),-np.ones(len(intervals),dtype=int)))
    order=np.argsort(points,kind='stable');points=points[order];active=np.cumsum(changes[order])
    widths=np.diff(points);area=float(np.dot(active[:-1],widths));elapsed=float(points[-1]-points[0])
    return dict(attempts=states,first_claim_utc=float(points[0]),last_accept_utc=float(points[-1]),
        elapsed_seconds=elapsed,peak_claimed_tasks=int(active.max()),
        average_claimed_tasks=area/elapsed,claimed_task_seconds=area,
        seconds_with_claimed_work=float(widths[active[:-1]>0].sum()),
        scope='Recorded claim-to-accept concurrency, including publication and dispatch. Recovery gaps are included in elapsed span. Not simultaneous CPU utilization.')


if __name__=='__main__':
    main(*sys.argv[1:])
