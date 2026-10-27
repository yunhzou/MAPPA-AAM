"""Derive seed-prefix and direction metrics exclusively from saved evidence."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import tarfile

from .evidence import read,save
from publication_analysis import union_outcome,representative_metrics,certificate_id


def archive_results(path):
    if not path.exists():return None
    with tarfile.open(path) as tar:
        manifest=json.load(tar.extractfile('archive-manifest.json'))['files']
        def checked(name):
            data=tar.extractfile(name).read()
            assert hashlib.sha256(data).hexdigest()==manifest[name]['sha256'],name
            return json.loads(data)
        result=dict(search=checked('search.json') if 'search.json' in manifest else {},seeds={})
        for name in manifest:
            if name.startswith('analysis/') and name.endswith('/evaluation.json'):
                row=checked(name);key=(row['cut'],row['seed'])
                assert key not in result['seeds'], ('duplicate evaluation', key)
                assert tuple(map(int,name.split('/')[1:3]))==(row['cut'],row['seed']-1)
                for field,member in [('raw_sha256',f'cuts/cut_{row["cut"]:05d}.raw.pkl.gz'),
                                     ('finalized_sha256',f'finalized/{row["cut"]:05d}.pkl.gz')]:
                    if field in row:assert row[field]==manifest[member]['sha256'], (name,field)
                result['seeds'][key]=dict(evaluation=row,classes=checked(name.replace('evaluation.json','classes.json')))
        return result


def merge(rows):
    groups={}
    for row in rows:
        origin=dict(direction=row['direction'],cut=row['cut'],seed=row['seed'],terminals=row['members'])
        old=groups.get(row['id'])
        if old is None:groups[row['id']]={**row,'origins':[origin]}
        else:
            origins=old['origins']+[origin]
            chosen=min((old,row),key=lambda x:(x['key'],x['direction'],x['cut'],x['seed']))
            groups[row['id']]={**chosen,'origins':origins}
    return sorted(groups.values(),key=lambda x:(x['key'],x['id']))


def case_metrics(payload):
    run,case,task_ids=payload
    run=Path(run);index=case['index'];modes={}
    if case['status']!='prepared':return dict(**case,modes={})
    reference_bytes=(run/'inputs'/str(index)/'reference.json').read_bytes()
    if 'reference_sha256' in case:
        assert hashlib.sha256(reference_bytes).hexdigest()==case['reference_sha256']
    reference=json.loads(reference_bytes)
    expected=certificate_id(**{k:v for k,v in reference.items() if k in ('features','mapping')})
    directed={d:archive_results(run/'archives'/f'{task_ids[d]:05d}.tar') for d in ('R_to_P','P_to_R')}
    smaller='R_to_P' if case['reactant_atoms']<=case['product_atoms'] else 'P_to_R'
    for scope in ('sweep','uncut'):
        for n in range(1,11):
            for mode,directions in [('R_to_P',['R_to_P']),('P_to_R',['P_to_R']),('smaller_first',[smaller]),('bidirectional',['R_to_P','P_to_R'])]:
                verdicts=[];classes=[];complete=True
                for d in directions:
                    data=directed[d]
                    count=data['search'].get('metrics',{}).get('cut_count') if data else None
                    selected=[r for (cut,seed),r in data['seeds'].items() if seed<=n and (scope=='sweep' or cut==0)] if data else []
                    keys={(cut,seed) for cut,seed in data['seeds'] if seed<=n and (scope=='sweep' or cut==0)} if data else set()
                    wanted=range(count) if scope=='sweep' and count is not None else range(1)
                    full=count is not None and keys=={(cut,seed) for cut in wanted for seed in range(1,n+1)}
                    outcomes=[r['evaluation']['reference_recovery'] for r in selected]
                    if not full:outcomes.append('unknown')
                    verdicts.append(union_outcome(outcomes));complete &= full
                    classes.extend(c for r in selected for c in r['classes'])
                key=f'{mode}/seeds{n}' if scope=='sweep' else f'uncut/{mode}/seeds{n}'
                ranked=merge(classes)
                modes[key]=dict(reference_recovery=union_outcome(verdicts),ranking_complete=complete,
                    representative=representative_metrics(ranked,expected),top10_classes=ranked[:10])
    coverage={d:dict(expected_cuts=data['search'].get('metrics',{}).get('cut_count'),
        saved_evaluation_rows=len(data['seeds']),unknown_rows=sum(r['evaluation']['reference_recovery']=='unknown' for r in data['seeds'].values()))
        if data else None for d,data in directed.items()}
    return dict(**case,modes=modes,coverage=coverage)


def summarize(run,output,workers=1):
    import concurrent.futures
    run,output=Path(run),Path(output);manifest=read(run/'manifest.json')
    tasks={(r['index'],r['direction']):n for n,r in enumerate(manifest['tasks'])}
    keys=[prefix+f'{mode}/seeds{n}' for prefix in ('','uncut/') for n in range(1,11)
          for mode in ('R_to_P','P_to_R','smaller_first','bidirectional')]
    totals={key:Counter() for key in keys}
    payloads=[(str(run),case,{d:tasks[case['index'],d] for d in ('R_to_P','P_to_R')}
               if case['status']=='prepared' else {}) for case in manifest['inputs']]
    def consume(records):
        for number,record in enumerate(records,1):
            for key,counter in totals.items():
                row=record['modes'].get(key)
                if row is None:counter['unknown']+=1;continue
                counter[row['reference_recovery']]+=1
                counter['ranking_complete']+=row['ranking_complete']
                if row['ranking_complete']:
                    metric=row['representative']
                    for k,hit in metric['topk'].items():counter['top'+k]+=hit
                    counter['top5_tie_inclusive']+=metric['top5_tie_inclusive']['recovered']
                    for gap,item in metric['event_windows'].items():counter['representative_window_'+gap]+=item['recovered']
            save(output/'cases'/f'{record["index"]:05d}.json',record)
            if number%100==0:print(f'{number}/{len(payloads)} cases summarized',flush=True)
    if workers==1:consume(map(case_metrics,payloads))
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            consume(pool.map(case_metrics,payloads,chunksize=1))
    summary=dict(total=manifest['total'],modes={k:dict(v) for k,v in totals.items()},
        scope='Raw compressed-family recovery and terminal-representative ranking; no claim of completed event decoding.',
        timing_scope=manifest['timing'],source_manifest_sha256=hashlib.sha256((run/'manifest.json').read_bytes()).hexdigest())
    save(output/'summary.json',summary)
    return summary


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('output');p.add_argument('--workers',type=int,default=1)
    a=p.parse_args();summarize(a.run,a.output,a.workers)
