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
    from mappa.event_patterns import SignedEventIndex
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
    normalized_rows=[];records=[];times=[]
    for item in payload['outputs']:
        value=read_blob(Blob(**item['blob']),workspace.scratch)
        key=f'{value["task_id"]:05d}/cut/{value["cut"]:05d}/seed/{value["seed"]-1:02d}'
        if key!=item['key'] or value['source']['input.json']['sha256']!=payload['input']['sha256'] or value['source']['reference.json']['sha256']!=payload['reference']['sha256']:
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
    return dict(index=payload['index'],modes=modes,records=records,timing=times,
        input=payload['input'],reference=payload['reference'],
        calculation=dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu),
        scientific_audit_accepted=False,
        reference_policy='Original Golden projected heavy-atom correspondence with stereo/charge labels. Hydrogen-unannotated references do not uniquely specify an all-atom event count.')
