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


def summarize(run,output):
    run,output=Path(run),Path(output);manifest=read(run/'manifest.json')
    tasks={(r['index'],r['direction']):n for n,r in enumerate(manifest['tasks'])}
    totals={f'{mode}/seeds{n}':Counter() for n in range(1,11) for mode in ('R_to_P','P_to_R','smaller_first','bidirectional')}
    for case in manifest['inputs']:
        index=case['index'];modes={}
        if case['status']!='prepared':
            for counter in totals.values():counter['unknown']+=1
            save(output/'cases'/f'{index:05d}.json',case);continue
        expected=certificate_id(**{k:v for k,v in read(run/'inputs'/str(index)/'reference.json').items() if k in ('features','mapping')})
        directed={d:archive_results(run/'archives'/f'{tasks[index,d]:05d}.tar') for d in ('R_to_P','P_to_R')}
        smaller='R_to_P' if case['reactant_atoms']<=case['product_atoms'] else 'P_to_R'
        for n in range(1,11):
            for mode,directions in [('R_to_P',['R_to_P']),('P_to_R',['P_to_R']),('smaller_first',[smaller]),('bidirectional',['R_to_P','P_to_R'])]:
                verdicts=[];classes=[];complete=True
                for d in directions:
                    data=directed[d]
                    count=data['search'].get('metrics',{}).get('cut_count') if data else None
                    selected=[r for (cut,seed),r in data['seeds'].items() if seed<=n] if data else []
                    full=count is not None and len(selected)==count*n
                    outcomes=[r['evaluation']['reference_recovery'] for r in selected]
                    if not full:outcomes.append('unknown')
                    verdicts.append(union_outcome(outcomes));complete &= full
                    classes.extend(c for r in selected for c in r['classes'])
                key=f'{mode}/seeds{n}';outcome=union_outcome(verdicts)
                ranked=merge(classes)
                metric=representative_metrics(ranked,expected)
                record=dict(reference_recovery=outcome,ranking_complete=complete,representative=metric,
                    top10_classes=ranked[:10])
                modes[key]=record;totals[key][outcome]+=1
                totals[key]['ranking_complete']+=complete
                if complete:
                    for k,hit in metric['topk'].items():totals[key]['top'+k]+=hit
        save(output/'cases'/f'{index:05d}.json',dict(**case,modes=modes))
    summary=dict(total=manifest['total'],modes={k:dict(v) for k,v in totals.items()},
        scope='Raw compressed-family recovery and terminal-representative ranking; no claim of completed event decoding.',
        timing_scope=manifest['timing'])
    save(output/'summary.json',summary)
    return summary


if __name__=='__main__':
    import sys
    summarize(sys.argv[1],sys.argv[2])
