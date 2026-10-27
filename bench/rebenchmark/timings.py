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


if __name__=='__main__':
    main(*sys.argv[1:])
