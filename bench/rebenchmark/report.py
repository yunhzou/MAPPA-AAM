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


def aggregate_fresh(plan, output, scratch):
    """Join verified case reports and compressed membership over the full grid."""
    from collections import defaultdict
    from concurrent.futures import ThreadPoolExecutor
    import gzip, sqlite3
    from mappa.mappa_core.core.scheduler import Blob,read_blob
    from mappa.mappa_core.utils.io.checkpoints import _loads_checkpoint
    from .timings import distribution,case_timing,execution_timeline
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    if sha(plan['reports'])!=plan['reports_sha256']:
        raise ValueError('Accepted report manifest changed')
    reports=read(plan['reports'])
    for key,item in plan.get('report_overrides',{}).items():
        if key in reports:raise ValueError('Completed report would be repeated')
        reports[key]=item
    expected={f'metrics/{i:05d}' for i in range(1851)}
    if set(reports)!=expected:
        raise ValueError('Final report lacks expected cases')
    if sha(plan['membership_queue'])!=plan['membership_queue_sha256']:
        raise ValueError('Membership checkpoint changed')
    db=sqlite3.connect('file:'+plan['membership_queue']+'?mode=ro&immutable=1',uri=True)
    members={k:(state,_loads_checkpoint(raw) if raw else None) for k,state,raw in db.execute('SELECT key,state,result FROM tasks')};db.close()
    if set(members)!={f'membership/{i:05d}' for i in range(1851)}:
        raise ValueError('Membership denominator differs')
    def extract(item):
        key,source=item;record=read_blob(Blob(**source['blob']),scratch)
        index=record['index']
        if key!=f'metrics/{index:05d}' or not record['scientific_audit_accepted'] or record['missing_tasks']:
            raise ValueError('Unverified or incomplete report inventory')
        state,resource=members[f'membership/{index:05d}']
        member=read_blob(resource['blob'],scratch) if state=='complete' else None
        if member is not None:
            if member['index']!=index or member['input']!=record['input'] or member['reference']!=record['reference']:
                raise ValueError('Membership belongs to another input')
            keys=sorted(r['key'] for r in record['records'])
            if len(keys)!=len(set(keys)) or hashlib.sha256(json.dumps(keys).encode()).hexdigest()!=member['expected_keys_sha256']:
                raise ValueError('Membership/report expected-key join differs')
        return record,member,case_timing(record),state,resource
    modes=defaultdict(Counter);phases=defaultdict(lambda:defaultdict(list))
    case_phases=defaultdict(lambda:defaultdict(list));direction_phases=defaultdict(list)
    prefix_times=defaultdict(Counter);counts=Counter();cases=[];membership_cpu=[];membership_wall=[];membership_peaks=[]
    query_phases=defaultdict(lambda:defaultdict(list))
    with ThreadPoolExecutor(max_workers=4) as workers:
        for record,member,timing,state,resource in workers.map(extract,sorted(reports.items())):
            counts.update(cases=1,conversions=len(record['conversion_checks']),search_outputs=record['search_outputs_verified'],
                incomplete_conversions=sum(not r['complete'] for r in record['conversion_checks']))
            for name,value in record['modes'].items():
                counter=modes[name];counter['cases']+=1
                counter['representative_recovered']+=value['representative_hits']['strict']
                counter['relaxed_representative_recovered']+=value['representative_hits']['without_both']
                counter['complete_windows']+=value['decoded_window_complete']
                for k,hit in value['representative_topk'].items():counter['representative_top'+k]+=hit
                for kind in ('aam','converted'):
                    status=member['modes'][name][kind]['status'] if member else 'unknown'
                    counter[kind+'_'+status]+=1
            for name,samples in timing['phases'].items():
                for measure,values in samples.items():
                    phases[name][measure].extend(values)
                    case_phases[name][measure].append(max(values,default=0) if measure=='peak' else sum(values))
            for direction,values in timing['directions'].items():
                for name,total in values.items():direction_phases[name].append(total)
            for name,values in timing['prefixes'].items():prefix_times[name].update(values)
            if member is not None:
                membership_cpu.append(member['timing']['cpu_seconds']);membership_wall.append(member['timing']['wall_seconds'])
                # Recovered published outputs retain operation timings but may
                # lack the supervisor's memory sample. Missing is not zero.
                peak = resource.get('peak_tree_rss_bytes')
                if peak is not None:
                    membership_peaks.append(peak)
                for rows in member['rows'].values():
                    for row in rows:
                        for query in row['queries']:
                            for phase,values in query.get('timing',{}).items():
                                for name,value in values.items():query_phases[phase][name].append(value)
            cases.append(dict(index=record['index'],modes={name:dict(representative=value['representative_hits']['strict'],
                complete=value['decoded_window_complete'],minimum=value['observed_minimum'],
                membership=member['modes'][name] if member else None) for name,value in record['modes'].items()},
                original_time=timing['directions'],membership_state=state,
                membership_time=member['timing'] if member else None))
            if len(cases)%100==0:print(f'{len(cases)}/1851 reports joined',flush=True)
    if (counts['conversions'],counts['search_outputs'])!=(1775520,1065312):
        raise ValueError('Final scientific output denominator differs')
    def statistics(groups):return {name:{measure:distribution(values) for measure,values in samples.items()} for name,samples in groups.items()}
    summary=dict(counts=dict(counts),modes={name:dict(value) for name,value in modes.items()},
        original_time=dict(per_task=statistics(phases),per_case=statistics(case_phases),
            per_direction={k:distribution(v) for k,v in direction_phases.items()},
            prefix_totals={k:dict(v) for k,v in prefix_times.items()}),
        membership_time=dict(cpu=distribution(membership_cpu),wall=distribution(membership_wall),
            peak=distribution(membership_peaks),peak_measured_cases=len(membership_peaks),
            peak_missing_cases=len(membership_cpu)-len(membership_peaks),queries=statistics(query_phases)),
        execution=execution_timeline(plan['scientific_queue']),
        scope='1851 Golden cases. Membership uses exact heavy-atom correspondence modulo endpoint chemistry in saved compressed choices. Converted membership is before event-window filtering; representative accuracy uses minimum+5 windows. Incomplete conversion flags remain explicit.',
        timing_scope='Original measured CPU and summed worker wall time. Matching tasks contain two seeds: odd-prefix matching totals include the paired next seed. Finalization is shared across all ten seeds. No standalone seed latency or hidden-bijection optimum is inferred.',
        sources=plan)
    save(output/'summary.json',summary)
    with gzip.open(output/'cases.json.gz','wt') as stream:json.dump(cases,stream,separators=(',',':'))
    save(output/'accepted-reports.json',reports)
    lines=['# Post-refactor Golden benchmark','',summary['scope'],'',
        '| Seeds | AAM recovered | AAM unresolved | Converted recovered | Converted unresolved | Representative hits |',
        '|---:|---:|---:|---:|---:|---:|']
    for n in range(1,11):
        mode=modes[f'sweep/bidirectional/seeds{n}']
        lines.append(f"| {n} | {mode['aam_recovered']}/1851 | {mode['aam_unknown']} | {mode['converted_recovered']}/1851 | {mode['converted_unknown']} | {mode['representative_recovered']}/1851 |")
    lines+=['','## Original computation time','',
        '| Stage | CPU total (h) | CPU median (s/task) | CPU p95 | CPU maximum | Worker wall total (h) | Peak memory (GiB) |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for name in ('matching_batch','finalization','conversion_decoding','conversion_task','conversion_publication'):
        values=summary['original_time']['per_task'][name];c=values['cpu'];w=values['wall'];p=values['peak']
        lines.append(f"| {name} | {c['total']/3600:.3f} | {c['median']:.6f} | {c['p95']:.6f} | {c['maximum']:.6f} | {w['total']/3600:.3f} | {p['maximum']/2**30:.3f} |")
    lines+=['',summary['timing_scope'],'',f"Incomplete conversion outputs: {counts['incomplete_conversions']}/{counts['conversions']}. These are not negative mapping results."]
    (output/'report.md').write_text('\n'.join(lines)+'\n')
    return summary


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('output');p.add_argument('--workers',type=int,default=1)
    p.add_argument('--replacements',type=Path)
    a=p.parse_args();summarize(a.run,a.output,a.workers,read(a.replacements) if a.replacements else None)
