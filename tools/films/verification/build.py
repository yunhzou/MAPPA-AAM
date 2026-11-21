"""Build the verification film from a recorded MAPPA replay, without re-searching."""
from pathlib import Path
import argparse,gzip,json,hashlib
import numpy as np


def display_coordinates(xyz):
    points=np.asarray(xyz,float);points-=points.mean(0)
    _,_,axes=np.linalg.svd(points,full_matrices=False)
    rotation=axes.T
    if np.linalg.det(rotation)<0:rotation[:,-1]*=-1
    return (points@rotation).tolist()


def build(run,template,library,output,preparation):
    output.mkdir(parents=True,exist_ok=True)
    raw=(run/'trajectory.json.gz').read_bytes();trace=json.loads(gzip.decompress(raw))
    report=json.loads((run/'verification.json').read_text())
    assert report['status']=='verified_connectivity'
    path=trace['runs'][0]['paths'][0];frames=[]
    for i,f in enumerate(path['frames']):
        if f['kind'] not in ('seed_start','commit','locked','terminal'):continue
        c=f['candidates'][f['preferred']] if f['candidates'] else None
        m={**f['locked'],**(c['witness'] if c else {})}
        if f['kind'] in ('seed_start','commit'):
            assert len(m)==len(frames)+1
        frames.append(dict(kind=f['kind'],mapping=m,active=f['active'],edge=f.get('edge'),
                           source_event_index=i,candidate_count=len(f['candidates'])))
    n=report['atoms']
    assert len(trace['runs'][0]['graph']['transitions']) == 1
    assert path['mapping'] == report['mapping']
    assert len(path['mapping']) == n == len(trace['input']['product']['elements'])
    ids=[path['mapping'][str(i)] for i in range(n)]
    r=np.asarray(trace['input']['reactant']['wbo']) >= .5
    p=np.asarray(trace['input']['product']['wbo'])[np.ix_(ids,ids)] >= .5
    assert np.array_equal(r,p)
    assert sum(f['kind']=='commit' for f in frames)==n-1
    assert frames[-1]['mapping']==path['mapping']
    inp=trace['input']
    for k in ('reactant','product'):
        e=inp[k];e['display_coordinates']=display_coordinates(e['coordinates'])
        e['edges']=[[i,j] for i in range(n) for j in range(i+1,n) if e['wbo'][i][j]>=.5]
        coords=np.asarray(e['coordinates']);display=np.asarray(e['display_coordinates'])
        assert np.allclose(np.linalg.norm(coords[:,None]-coords[None,:],axis=2),np.linalg.norm(display[:,None]-display[None,:],axis=2))
    for f in frames:
        m={int(a):b for a,b in f['mapping'].items()}
        assert len(m)==len(set(m.values()))
        assert all(inp['reactant']['elements'][a]==inp['product']['elements'][b] for a,b in m.items())
    data=dict(title='Verify what your model built.',input=inp,frames=frames,
              report=report,preparation=json.loads(preparation.read_text()),
              trace_sha256=hashlib.sha256(raw).hexdigest(),duration=28)
    (output/'film-data.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    html=template.read_text().replace('__LIBRARY__',library.read_text()).replace('__DATA__',json.dumps(data,separators=(',',':')).replace('</','<\\/'))
    (output/'index.html').write_text(html)
    (output/'science-validation.json').write_text(json.dumps(dict(status='passed',atoms=n,fragments=1,
        unit_growth_events=n-1,seed_placements=1,missing_connections=0,extra_connections=0,
        trace_sha256=data['trace_sha256'],all_display_frames_are_recorded=True,
        display_transforms_preserve_distances=True,scope=report['scope']),indent=2)+'\n')
    print(output/'index.html')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--preparation',type=Path,required=True)
    args=parser.parse_args();build(args.run,Path(__file__).with_name('film.html'),args.library,args.output,args.preparation)
