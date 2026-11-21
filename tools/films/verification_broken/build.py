"""Build a truthful broken-molecule film from the selected diagnostic replay."""
from pathlib import Path
import argparse,gzip,hashlib,json
import numpy as np


def build(run,library,preparation,intact,output):
    output.mkdir(parents=True,exist_ok=True)
    trace_bytes=(run/'trajectory.json.gz').read_bytes();trace=json.loads(gzip.decompress(trace_bytes))
    report=json.loads((run/'verification.json').read_text());prep=json.loads(preparation.read_text())
    assert report['status']=='different_connectivity' and report['fragments']==2
    assert report['missing']==1 and report['extra']==0 and report['complete']
    path=trace['runs'][0]['paths'][0];inp=trace['input'];n=report['atoms'];m={int(i):j for i,j in path['mapping'].items()}
    assert path['mapping']==report['mapping'] and len(m)==len(set(m.values()))==n
    wr=np.asarray(inp['reactant']['wbo']);wp=np.asarray(inp['product']['wbo']);ids=[m[i] for i in range(n)]
    assert np.triu((wr>=.5)&~(wp[np.ix_(ids,ids)]>=.5),1).sum()==1
    assert np.triu(~(wr>=.5)&(wp[np.ix_(ids,ids)]>=.5),1).sum()==0
    events=path['events'];assert len(events)==1 and events[0]['kind']=='broken'
    assert sorted(events[0]['p'])==sorted(prep['candidate_edge'])
    frames=[];clock=4.5;groups={};growth=0
    for i,f in enumerate(path['frames']):
        if f['kind'] not in ('seed_start','commit','consumed','locked','terminal'):continue
        c=f['candidates'][f['preferred']] if f['candidates'] else None
        mapping={**f['locked'],**(c['witness'] if c else {})}
        if f['kind'] in ('seed_start','commit'):
            growth+=1;assert len(mapping)==growth
            for atom in f['active']:groups[atom]=f['stage']-1
        length=.085
        if f['kind']=='consumed':length=1.6 if f['stage']==1 else .4
        if f['kind']=='locked':length=.8 if f['stage']==1 else .3
        if f['kind']=='terminal':length=0
        frames.append(dict(kind=f['kind'],stage=f['stage'],mapping=mapping,active=f['active'],
            edge=f.get('edge'),source_event_index=i,start=round(clock,6),duration=length,growth_step=growth-1))
        clock+=length
    assert growth==135 and len(groups)==135
    assert sorted([list(groups.values()).count(i) for i in range(2)])==[21,114]
    # Each endpoint is rigidly posed for display. During the introduction only,
    # candidate atoms interpolate from the intact control to the fixed broken XYZ.
    intact_rows=intact.read_text().splitlines()[2:];old=np.array([[float(x) for x in line.split()[1:4]] for line in intact_rows])
    for side in ['reactant','product']:
        e=inp[side];coords=np.asarray(e['coordinates']);center=coords.mean(0);_,_,vt=np.linalg.svd(coords-center,full_matrices=False);rot=vt.T
        if np.linalg.det(rot)<0:rot[:,-1]*=-1
        e['display_coordinates']=((coords-center)@rot).tolist()
        e['edges']=[[a,b] for a in range(n) for b in range(a+1,n) if e['wbo'][a][b]>=.5]
        if side=='product':e['intact_display_coordinates']=((old-center)@rot).tolist()
        d=np.asarray(e['display_coordinates']);assert np.allclose(np.linalg.norm(coords[:,None]-coords[None,:],axis=2),np.linalg.norm(d[:,None]-d[None,:],axis=2))
    data=dict(title='Break one connection.',input=inp,frames=frames,groups=groups,report=report,preparation=prep,
              missing_edge=events[0],growth_end=round(clock,6),duration=30,trace_sha256=hashlib.sha256(trace_bytes).hexdigest())
    html=Path(__file__).with_name('film.html').read_text().replace('__LIBRARY__',library.read_text()).replace('__DATA__',json.dumps(data,separators=(',',':')).replace('</','<\\/'))
    (output/'index.html').write_text(html);(output/'film-data.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    (output/'science-validation.json').write_text(json.dumps(dict(status='passed',atoms=135,fragments=2,fragment_sizes=[114,21],
        seed_placements=2,unit_growth_events=133,recorded_boundary_deferrals=2,missing_connections=1,extra_connections=0,
        candidate_components=2,all_growth_frames_are_recorded=True,construction_animation_is_illustrative=True,
        actual_broken_candidate_edge=events[0]['p'],trace_sha256=data['trace_sha256']),indent=2)+'\n')
    print('Growth ends:',data['growth_end'],'; output:',output/'index.html')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['run','library','preparation','intact','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();build(a.run,a.library,a.preparation,a.intact,a.output)
