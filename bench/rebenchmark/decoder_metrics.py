"""Summarize saved decoded witnesses; proof-audit acceptance is a separate gate."""
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time


def feature_variants(features):
    """Keep the Golden projection policy; relax only the named scoring labels."""
    variants = {}
    for name, stereo, charge in [('strict',False,False), ('without_stereo',True,False),
                                  ('without_charge',False,True), ('without_both',True,True)]:
        value = deepcopy(features)
        for side in value:
            side['colors'] = [list(c) for c in side['colors']]
            for color in side['colors']:
                if stereo: color[4] = ''
                if charge: color[1] = 0
            if stereo:
                side['bonds'] = [(a,b,(color[0],'')) for a,b,color in side['bonds']]
        variants[name] = value
    return variants


def summarize(rows, expected, reference_ids, slack=5):
    """The full expected key set gates tier, minimum and window completeness."""
    expected = set(expected)
    if len({r['key'] for r in rows}) != len(rows) or any(r['key'] not in expected for r in rows):
        raise ValueError('Duplicate or unexpected decoder output')
    present = {r['key'] for r in rows}
    tier = max((tuple(r['coverage']) for r in rows if r['coverage'] is not None), default=None)
    selected = [r for r in rows if r['coverage'] is not None and tuple(r['coverage']) == tier]
    minimum = min((c['total'] for r in selected for c in r['candidates']), default=None)
    window = None if minimum is None else minimum+slack
    def covers(row, bound):
        return row['bound'] is None or (bound is not None and row['bound'] >= bound)
    full = present == expected
    tier_known = full and all(r['catalogue_complete'] for r in rows)
    complete_paths = bool(selected) and all(r['decoded_complete'] for r in selected)
    proven = tier_known and minimum is not None and (minimum == 0 or
        (complete_paths and all(covers(r,minimum) for r in selected)))
    complete = tier_known and minimum is not None and complete_paths and all(covers(r,window) for r in selected)
    events, classes = {}, {}
    for row in selected:
        for c in row['candidates']:
            if window is not None and c['total'] <= window:
                if c['event_id'] in events and events[c['event_id']] != c['total']:
                    raise ValueError('Canonical event identity has inconsistent counts')
                events[c['event_id']] = c['total']
                old = classes.get(c['strict'])
                if old is None or (c['total'],c['event_id']) < (old['total'],old['event_id']):
                    classes[c['strict']] = c
    ranked = sorted(classes.values(),key=lambda c:(c['total'],c['strict']))
    # Relaxed hits inspect all witnesses, not only a chosen strict-class witness.
    witnesses = [c for r in selected for c in r['candidates'] if window is not None and c['total'] <= window]
    reference_events = [c['total'] for c in witnesses if c['strict'] == reference_ids['strict']]
    tie_bound = ranked[min(5,len(ranked))-1]['total'] if ranked else None
    windows = {}
    for gap in range(slack+1):
        bound = None if minimum is None else minimum+gap
        eligible = [c for c in witnesses if bound is not None and c['total'] <= bound]
        windows[str(gap)] = dict(event_classes=sum(v <= bound for v in events.values()) if bound is not None else 0,
            strict_representative_hit=any(c['strict']==reference_ids['strict'] for c in eligible))
    return dict(expected_rows=len(expected),present_rows=len(rows),missing_rows=sorted(expected-present),
        coverage=tier,coverage_tier_known=tier_known,selected_rows=len(selected),observed_minimum=minimum,
        minimum_proven=proven,decoded_window_complete=complete,event_window=window,
        event_class_count=len(events),strict_representative_classes=len(ranked),windows=windows,
        representative_hits={name:any(c[name]==cid for c in witnesses) for name,cid in reference_ids.items()},
        representative_topk={str(k):any(c['strict']==reference_ids['strict'] for c in ranked[:k]) for k in (1,3,5,10)},
        top5_tie_inclusive=any(c['strict']==reference_ids['strict'] for c in ranked if c['total']<=tie_bound),
        reference_compatible_witness_minimum=min(reference_events,default=None),
        reference_compatible_witness_offset=None if not reference_events else min(reference_events)-minimum,
        incomplete_rows=[r['key'] for r in selected if not r['decoded_complete']],
        insufficient_window_rows=[r['key'] for r in selected if not covers(r,window)],
        scientific_audit_accepted=False,
        scope='Saved event-class representatives. Reference-compatible witness offset is not a uniquely annotated all-H ground-truth event count. Final acceptance requires matching proof-audit output identities.')


def case_metrics(payload, workspace):
    """Read accepted output blobs once and normalize concrete reverse witnesses."""
    from mappa import AAMProblem, MolecularEndpoint
    from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
    from mappa.mappa_core.core.scheduler import Blob, read_blob
    from publication_analysis import certificate_id
    started,cpu=time.perf_counter(),time.process_time()
    def read_input(name):
        item=payload[name];data=Path(item['path']).read_bytes()
        if hashlib.sha256(data).hexdigest()!=item['sha256']:
            raise ValueError('Input checksum differs')
        return json.loads(data)
    raw, reference = read_input('input'), read_input('reference')
    problem=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw['name'])
    variants=feature_variants(reference['features'])
    expected_ids={name:certificate_id(features,reference['mapping']) for name,features in variants.items()}
    event_index=SignedEventIndex(problem,threshold=.5,metal_threshold=.3)
    @lru_cache(maxsize=65536)
    def normalized(pairs):
        mapping=dict(pairs)
        event=event_index.describe(event_index.vector(mapping))
        return dict(event_id=event['id'],total=event['total'])
    @lru_cache(maxsize=65536)
    def certificates(heavy):
        # certificate_id takes original atom indices, so retain those here.
        return {name:certificate_id(features,dict(heavy)) for name,features in variants.items()}
    fresh = payload.get('source_kind') == 'fresh_conversion'
    verified=[]
    if fresh:
        from .check_conversion_outputs import verify_conversion
        from mappa.mappa_core.utils.io.records import problem_from_record
        @lru_cache(maxsize=2)
        def fresh_input(blob):
            return problem_from_record(read_blob(blob,workspace.scratch))
        @lru_cache(maxsize=2)
        def finalized(blob):
            return read_blob(blob,workspace.scratch)['graph']
        for item in payload['search_outputs']:
            saved=read_blob(Blob(**item['blob']),workspace.scratch)
            if not saved['graph'].contexts:
                raise ValueError('Saved search output lacks contexts')
    normalized_rows=[];records=[];times=[]
    for item in payload['outputs']:
        value=read_blob(Blob(**item['blob']),workspace.scratch)
        key=f'{value["task_id"]:05d}/cut/{value["cut"]:05d}/seed/{value["seed"]-1:02d}'
        if key!=item['key']:
            raise ValueError('Decoded output task identity differs')
        if fresh:
            expected=item['payload']
            if (value['input']!=expected['input'] or value['finalized_sha256']!=value['finalized'].sha256
                    or value['finalized']!=Blob(**item['finalized_blob'])):
                raise ValueError('Fresh conversion source differs')
            if (value['cut'],value['seed'])!=(expected['number'],expected['seed']+1):
                raise ValueError('Fresh conversion seed/cut differs')
            directed=fresh_input(value['input'])
            # Verify the directed saved endpoint against the hashed original input.
            original=(problem.product,problem.reactant) if payload['directions'][str(value['task_id'])]=='P_to_R' else (problem.reactant,problem.product)
            import numpy as np
            for actual,wanted in zip((directed.reactant,directed.product),original):
                if tuple(actual.elements)!=tuple(wanted.elements) or not np.array_equal(actual.wbo,wanted.wbo):
                    raise ValueError('Fresh endpoint differs from original input')
            if not finalized(value['finalized']).contexts:
                raise ValueError('Finalized source lacks contexts')
            verified.append(verify_conversion(directed,value,item['resource'],key))
        elif value['source']['input.json']['sha256']!=payload['input']['sha256'] or value['source']['reference.json']['sha256']!=payload['reference']['sha256']:
            raise ValueError('Decoded output source or task identity differs')
        config=value['decoded']['config']
        if (config['threshold'],config['metal_threshold'])!=(.5,.3):
            raise ValueError('Unexpected event thresholds')
        direction=payload['directions'][str(value['task_id'])]
        candidates=[]
        for c in value['decoded']['candidates']:
            mapping={int(a):int(b) for a,b in c['mapping'].items()}
            if direction=='P_to_R':mapping={b:a for a,b in mapping.items()}
            pairs=tuple(sorted(mapping.items()));event=normalized(pairs)
            if event['total']!=c['total']:
                raise ValueError('Direction-normalized bond count differs')
            heavy=tuple((a,b) for a,b in pairs if problem.reactant.elements[a]!='H' and problem.product.elements[b]!='H')
            candidates.append(dict(event,**certificates(heavy)))
        decoded=value['decoded']
        normalized_rows.append(dict(key=key,direction=direction,cut=value['cut'],seed=value['seed'],
            coverage=value['coverage'],candidates=candidates,decoded_complete=decoded['complete'],
            catalogue_complete=decoded['incomplete_paths']==0,
            bound=decoded['event_window']))
        records.append(dict(key=key,blob=item['blob'],decoded_complete=decoded['complete'],
            incomplete_families=[r['family_id'] for r in decoded['family_reports'] if not r['complete']]))
        times.append(dict(key=key,**value['timing']))
    modes={}
    smaller='R_to_P' if problem.source_atom_count<=problem.target_atom_count else 'P_to_R'
    for scope in ('uncut','sweep'):
        for n in range(1,11):
            for mode,directions in [('R_to_P',{'R_to_P'}),('P_to_R',{'P_to_R'}),
                                    ('smaller_first',{smaller}),('bidirectional',{'R_to_P','P_to_R'})]:
                def included(key):
                    task,_,cut,_,seed=key.split('/')
                    return payload['directions'][str(int(task))] in directions and int(seed)<n and (scope=='sweep' or int(cut)==0)
                expected=[k for k in payload['expected'] if included(k)]
                rows=[r for r in normalized_rows if included(r['key'])]
                modes[f'{scope}/{mode}/seeds{n}']=summarize(rows,expected,expected_ids)
    if fresh:
        for mode in modes.values():
            mode['scientific_audit_accepted']=True
            mode['scope']='Independently counted retained converted branches; complete flags and missing keys preserved. No hidden-bijection or global-minimum claim.'
    return dict(index=payload['index'],modes=modes,records=records,timing=times,
        input=payload['input'],reference=payload['reference'],
        calculation=dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu),
        scientific_audit_accepted=fresh,conversion_checks=verified,
        search_outputs_verified=len(payload.get('search_outputs',[])),
        missing_tasks=payload.get('missing_tasks',[]),
        search_resources=[dict(key=r['key'],**r['resource']) for r in payload.get('search_outputs',[])],
        reference_policy='Original Golden projected heavy-atom correspondence with stereo/charge labels. Hydrogen-unannotated references do not uniquely specify an all-atom event count.')


def membership_summary(rows, expected):
    """An existential hit survives missing work; a negative needs every key."""
    expected = set(expected)
    selected = [r for r in rows if r['key'] in expected]
    hits = [r for r in selected if r['status'] == 'recovered']
    if hits:
        first = min(hits, key=lambda r: (r['seed'], r['cut'], r['key']))
        return dict(status='recovered', witness_key=first['key'], expected_rows=len(expected))
    present = {r['key'] for r in selected}
    status = 'not_recovered' if present == expected and all(r['status']=='not_recovered' for r in selected) else 'unknown'
    return dict(status=status, expected_rows=len(expected), checked_rows=len(selected),
                unchecked_rows=len(expected-present))


def case_membership(payload, workspace):
    """Check saved AAM and converted choices, stopping each scope at a witness.

    Converted coverage concerns the complete retained catalogue, before event
    window filtering. Windowed representative accuracy is reported separately.
    No matching, conversion or bijection enumeration is performed here.
    """
    from mappa import AAMProblem
    from mappa.mappa_core.core.scheduler import Blob, read_blob
    from mappa.mappa_core.utils.io.records import problem_from_record
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
    from golden_evaluation import reference_membership
    from .search import select_contexts
    start, cpu = time.perf_counter(), time.process_time()
    inputs={}
    for name in ('input','reference'):
        item=payload[name]; data=Path(item['path']).read_bytes()
        if hashlib.sha256(data).hexdigest()!=item['sha256']:
            raise ValueError('Membership input checksum differs')
        inputs[name]=json.loads(data)
    original=problem_from_record(inputs['input'])
    reference={int(a):int(b) for a,b in dict(inputs['reference']['mapping']).items()}
    features=inputs['reference']['features']
    outputs={item['key']:item for item in payload['outputs']}
    finals={item['key']:item for item in payload['search_outputs'] if item['key'].endswith('/finalize')}
    rows={'aam':[], 'converted':[]}; reads={'aam':0,'converted':0}
    @lru_cache(maxsize=1)
    def finalized(blob):
        reads['aam']+=1
        return read_blob(blob,workspace.scratch)['graph']
    for direction, label in payload['directions'].items():
        reverse=label=='P_to_R'
        problem=AAMProblem(original.product,original.reactant,original.name) if reverse else original
        ordered=list(reversed(features)) if reverse else features
        target={b:a for a,b in reference.items()} if reverse else reference
        keys=sorted((k for k in payload['expected'] if int(k.split('/')[0])==int(direction)),
                    key=lambda k:(int(k.split('/')[-1]),int(k.split('/')[2])))
        for kind in rows:
            deadline=time.perf_counter()+payload.get('membership_seconds',120.)
            found=set()
            for key in keys:
                _,_,cut,_,seed=key.split('/');cut,seed=int(cut),int(seed)+1
                scopes={'sweep','uncut'} if cut==0 else {'sweep'}
                if scopes<=found:
                    continue
                if time.perf_counter()>=deadline:
                    break
                if kind=='aam':
                    parent=key.rsplit('/seed/',1)[0]+'/finalize'
                    item=finals.get(parent)
                    if item is None:
                        continue
                    blob=Blob(**item['blob'])
                    graph=select_contexts(finalized(blob),[seed-1])
                    representatives=((t,graph.states[t].mapping) for t in graph.terminals)
                    paths=lambda:((list(path.transitions),path) for path in graph.paths())
                    complete=True  # membership is scoped to this retained graph, including capped searches
                else:
                    item=outputs.get(key)
                    if item is None:
                        continue
                    blob=Blob(**item['blob']);value=read_blob(blob,workspace.scratch);reads[kind]+=1
                    if (value['task_id'],value['cut'],value['seed'])!=(int(direction),cut,seed):
                        raise ValueError('Converted membership task identity differs')
                    catalogue=FinalBranchCatalogue.from_record(problem,value['decoded']['converted_catalogue'])
                    representatives=((i,f.mapping) for i,f in enumerate(catalogue.families))
                    paths=lambda:((i,f.as_path(problem)) for i,f in enumerate(catalogue.families))
                    complete=value['decoded']['complete']
                answer=reference_membership(problem,representatives,paths,ordered,target,
                                             seconds=max(0,deadline-time.perf_counter()))
                if answer['witness'] is not None:
                    mapping=dict(answer['witness']['mapping'])
                    answer['witness']['input_orientation_mapping']=sorted(
                        ({b:a for a,b in mapping.items()} if reverse else mapping).items())
                rows[kind].append(dict(key=key,seed=seed,cut=cut,direction=label,source_blob=item['blob'],
                                       conversion_complete=complete,**answer))
                if answer['status']=='recovered':
                    found.update(scopes)
    modes={};smaller='R_to_P' if original.source_atom_count<=original.target_atom_count else 'P_to_R'
    for scope in ('uncut','sweep'):
        for n in range(1,11):
            for mode,directions in [('R_to_P',{'R_to_P'}),('P_to_R',{'P_to_R'}),
                                    ('smaller_first',{smaller}),('bidirectional',{'R_to_P','P_to_R'})]:
                expected=[key for key in payload['expected']
                    if payload['directions'][str(int(key.split('/')[0]))] in directions
                    and int(key.split('/')[-1])<n and (scope=='sweep' or int(key.split('/')[2])==0)]
                modes[f'{scope}/{mode}/seeds{n}']={kind:membership_summary(value,expected) for kind,value in rows.items()}
    return dict(index=payload['index'],modes=modes,rows=rows,reads=reads,
        input=payload['input'],reference=payload['reference'],
        expected_keys_sha256=hashlib.sha256(json.dumps(sorted(payload['expected'])).encode()).hexdigest(),
        expected_keys=len(payload['expected']),missing_tasks=payload['missing_tasks'],
        timing=dict(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu),
        scope='Exact Golden heavy-atom relation modulo endpoint chemical symmetry in saved AAM and all retained converted choices, before event-window filtering. No global search-completeness claim.')
