"""Derive seed-prefix and direction metrics exclusively from saved evidence."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import tarfile

from .evidence import read,save,sha
from publication_analysis import union_outcome,representative_metrics,certificate_id


def archive_results(path, overlays=()):
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
        for saved in overlays:
            row = saved['evaluation']
            key = (row['cut'], row['seed'])
            assert saved['state'] == 'complete' and key not in result['seeds'], key
            for name, descriptor in saved['source'].items():
                assert descriptor['sha256'] == manifest[name]['sha256'], name
            raw = f'cuts/cut_{row["cut"]:05d}.raw.pkl.gz'
            assert row['raw_sha256'] == saved['raw']['sha256'] == manifest[raw]['sha256']
            final = f'finalized/{row["cut"]:05d}.pkl.gz'
            if final in manifest:
                assert row['finalized_sha256'] == manifest[final]['sha256']
            else:
                assert saved['finalized'] is not None
                assert row['finalized_sha256'] == saved['finalized']['sha256']
            previous_classes = f'analysis/{row["cut"]:05d}/{row["seed"]-1:02d}/classes.json'
            if previous_classes in manifest:
                assert saved['classes'] == checked(previous_classes), previous_classes
            result['seeds'][key] = dict(evaluation=row, classes=saved['classes'])
        return result


def load_overlays(specification):
    """Accept hashed continuation records once, with no duplicate task/seed keys."""
    outputs, seen = {}, set()
    for source in specification:
        assert sha(source['path']) == source['sha256']
        for key, record in read(source['path']).items():
            if record['state'] != 'complete':
                continue
            task, _, cut, _, seed = key.split('/')
            row = record['evaluation']
            identity = (int(task), int(cut), int(seed)+1)
            assert identity not in seen
            assert identity[1:] == (row['cut'], row['seed'])
            seen.add(identity)
            outputs.setdefault(int(task), []).append(record)
    return outputs


def resolve_archives(run, manifest, replacements=None):
    """Accept a continuation only for the same task with all saved cuts intact."""
    paths={n:run/'archives'/f'{n:05d}.tar' for n in range(len(manifest['tasks']))}
    if replacements is None:
        return paths
    assert replacements['source_manifest_sha256']==sha(run/'manifest.json')
    seen=set()
    for row in replacements['replacements']:
        task=row['task']
        assert task not in seen and task in paths
        seen.add(task)
        replacement=Path(row['archive'])
        assert sha(replacement)==row['sha256']
        with tarfile.open(paths[task]) as old, tarfile.open(replacement) as new:
            before=json.load(old.extractfile('archive-manifest.json'))['files']
            after=json.load(new.extractfile('archive-manifest.json'))['files']
            data=new.extractfile('task.json').read()
            assert hashlib.sha256(data).hexdigest()==after['task.json']['sha256']
            recorded=json.loads(data)
            for field in ('index','direction','config'):
                assert recorded[field]==manifest['tasks'][task][field], field
            for name,entry in before.items():
                if name.startswith('cuts/') and name.endswith('.raw.pkl.gz'):
                    assert after[name]==entry, name
        paths[task]=replacement
    return paths


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
    run,case,task_ids,*extra=payload
    overlays = extra[0] if extra else {}
    run=Path(run);index=case['index'];modes={}
    if case['status']!='prepared':return dict(**case,modes={})
    reference_bytes=(run/'inputs'/str(index)/'reference.json').read_bytes()
    if 'reference_sha256' in case:
        assert hashlib.sha256(reference_bytes).hexdigest()==case['reference_sha256']
    reference=json.loads(reference_bytes)
    expected=certificate_id(**{k:v for k,v in reference.items() if k in ('features','mapping')})
    directed={d:archive_results(Path(task_ids[d]), overlays.get(d, ())) for d in ('R_to_P','P_to_R')}
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


def count_record(record, counters, sign=1):
    for key, counter in counters.items():
        row = record['modes'].get(key)
        if row is None:
            counter['unknown'] += sign
            continue
        counter[row['reference_recovery']] += sign
        counter['ranking_complete'] += sign * row['ranking_complete']
        if row['ranking_complete']:
            metric = row['representative']
            for k, hit in metric['topk'].items():
                counter['top'+k] += sign * hit
            counter['top5_tie_inclusive'] += sign * metric['top5_tie_inclusive']['recovered']
            for gap, item in metric['event_windows'].items():
                counter['representative_window_'+gap] += sign * item['recovered']


def refresh(source, run, output, specifications):
    """Update only cases affected by validated new checks; preserve the base report."""
    run, output = Path(run), Path(output)
    manifest = read(run/'manifest.json')
    overlays = load_overlays(specifications)
    changed = {manifest['tasks'][task]['index'] for task in overlays}
    task_ids = {(row['index'], row['direction']): task for task, row in enumerate(manifest['tasks'])}
    with tarfile.open(source) as tar:
        members = json.load(tar.extractfile('archive-manifest.json'))['files']
        def checked(name):
            data = tar.extractfile(name).read()
            assert hashlib.sha256(data).hexdigest() == members[name]['sha256']
            return json.loads(data)
        summary = checked('summary.json')
        assert summary['source_manifest_sha256'] == sha(run/'manifest.json')
        counters = {k: Counter(v) for k, v in summary['modes'].items()}
        archives = resolve_archives(run, manifest, summary.get('replacements'))
        for case in manifest['inputs']:
            index = case['index']
            if index not in changed:
                continue
            ids = {d: task_ids[index, d] for d in ('R_to_P', 'P_to_R')}
            paths = {d: str(archives[task]) for d, task in ids.items()}
            additions = {d: overlays.get(task, ()) for d, task in ids.items()}
            before = checked(f'cases/{index:05d}.json')
            after = case_metrics((str(run), case, paths, additions))
            count_record(before, counters, -1)
            count_record(after, counters)
            save(output/'cases'/f'{index:05d}.json', after)
        summary.update(modes={k:dict(v) for k,v in counters.items()},
            continued_checks=specifications, changed_cases=sorted(changed),
            base_report=dict(path=str(source), sha256=sha(source)))
        assert all(v >= 0 for counts in counters.values() for v in counts.values())
        save(output/'summary.json', summary)
        return summary


def summarize(run,output,workers=1,replacements=None):
    import concurrent.futures
    run,output=Path(run),Path(output);manifest=read(run/'manifest.json')
    archives=resolve_archives(run,manifest,replacements)
    tasks={(r['index'],r['direction']):str(archives[n]) for n,r in enumerate(manifest['tasks'])}
    keys=[prefix+f'{mode}/seeds{n}' for prefix in ('','uncut/') for n in range(1,11)
          for mode in ('R_to_P','P_to_R','smaller_first','bidirectional')]
    totals={key:Counter() for key in keys}
    payloads=[(str(run),case,{d:tasks[case['index'],d] for d in ('R_to_P','P_to_R')}
               if case['status']=='prepared' else {}) for case in manifest['inputs']]
    def consume(records):
        for number,record in enumerate(records,1):
            count_record(record, totals)
            save(output/'cases'/f'{record["index"]:05d}.json',record)
            if number%100==0:print(f'{number}/{len(payloads)} cases summarized',flush=True)
    if workers==1:consume(map(case_metrics,payloads))
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            consume(pool.map(case_metrics,payloads,chunksize=1))
    summary=dict(total=manifest['total'],modes={k:dict(v) for k,v in totals.items()},
        scope='Raw compressed-family recovery and terminal-representative ranking; no claim of completed event decoding.',
        replacements=replacements, timing_scope=manifest['timing'],source_manifest_sha256=hashlib.sha256((run/'manifest.json').read_bytes()).hexdigest())
    save(output/'summary.json',summary)
    return summary


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('output');p.add_argument('--workers',type=int,default=1)
    p.add_argument('--replacements',type=Path)
    a=p.parse_args();summarize(a.run,a.output,a.workers,read(a.replacements) if a.replacements else None)
