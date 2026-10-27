"""Independently rescore saved witnesses against original NPZ pairs."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np


def rescore(elements,r,p,mapping):
    # Ir is the only metal present in either supplied example.
    assert set(elements)<= {'C','H','I','Ir','O','S','F','N'}
    n=len(elements); pairs=np.triu_indices(n,1)
    vector=np.array([mapping[i] for i in range(n)])
    delta=p[vector[pairs[0]],vector[pairs[1]]]-r[pairs]
    threshold=np.array([.3 if 'Ir' in (elements[i],elements[j]) else .5 for i,j in zip(*pairs)])
    return {kind:sorted((int(i),int(j)) for i,j in zip(pairs[0][mask],pairs[1][mask])) for kind,mask in [('broken',delta<=-threshold),('formed',delta>=threshold)]}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--result',type=Path,required=True);parser.add_argument('--inputs',type=Path,required=True);args=parser.parse_args()
    root=args.result;cfg=json.loads((root/'configuration.json').read_text());search=json.loads((root/'search.json').read_text())
    data=json.load(gzip.open(root/'decoded.json.gz','rt'));catalogue=json.load(gzip.open(root/'catalogue.json.gz','rt'))
    arrays=[]
    for name in ['reference_arrays.npz','prediction_arrays.npz']:
        path=args.inputs/name
        assert hashlib.sha256(path.read_bytes()).hexdigest()==cfg['inputs'][name]
        with np.load(path,allow_pickle=False) as npz:arrays.append({k:npz[k].copy() for k in npz.files})
    r,p=arrays;n=len(r['elements']); hist={}
    for candidate in data['candidates']:
        mapping={int(a):int(b) for a,b in candidate['mapping'].items()}
        assert sorted(mapping)==list(range(n)) and sorted(mapping.values())==list(range(n))
        assert all(r['elements'][a]==p['elements'][b] for a,b in mapping.items())
        events=rescore(r['elements'],r['wbo'],p['wbo'],mapping)
        assert events=={k:sorted(tuple(x) for x in v) for k,v in candidate['events'].items()}
        assert candidate['total']==sum(map(len,events.values()))
        hist[candidate['total']]=hist.get(candidate['total'],0)+1
    reports=data['family_reports']
    assert [r['family_id'] for r in reports]==list(range(len(catalogue['families'])))
    complete=all(r['complete'] for r in reports) and data['incomplete_paths']==0
    assert data['complete']==complete
    assert data['config']['max_events'] is None and data['config']['event_slack'] is None
    output=dict(case=cfg['case'],atoms=n,complete=complete,cuts=search['cuts'],search_capped=search['capped'],family_count=len(reports),complete_family_count=sum(r['complete'] for r in reports),incomplete_paths=data['incomplete_paths'],candidate_count=len(data['candidates']),event_histogram=hist,minimum_events=min(hist,default=None),all_witnesses_independently_rescored=True,input_hashes_verified=True,unrestricted_event_window=True,geometry_selection=False,evidence_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ['decoded.json.gz','catalogue.json.gz','search.json','configuration.json']})
    (root/'independent-audit.json').write_text(json.dumps(output,indent=2));print(json.dumps(output,indent=2))

if __name__=='__main__':main()
