"""Reference-blind class ranking and separate publication accuracy queries."""
import hashlib

import pynauty
from golden_evaluation import colored_graph, project


def certificate_id(features, mapping):
    return hashlib.sha256(pynauty.certificate(colored_graph(features,project(mapping,features)))).hexdigest()


def rank_archive(aam, plan):
    """No reference input. Canonical classes in original R/P space."""
    features = [dict(e.metadata)['scoring_features'] for e in
                (plan.input_problem.reactant,plan.input_problem.product)]
    score = ranker(plan.input_problem)
    cache, groups = {}, {}
    for terminal in aam.graph.terminals:
        mapping = plan.to_input_mapping(aam.graph.states[terminal].mapping)
        heavy = tuple(sorted(project(mapping,features).items()))
        if heavy not in cache:
            cache[heavy] = certificate_id(features,mapping)
        cid = cache[heavy]
        key, events = score(mapping)
        if cid not in groups:
            groups[cid] = dict(id=cid,key=key,events=events,terminal=terminal,
                              direction=plan.direction,members=[])
        row = groups[cid]
        row['members'].append(terminal)
        if key < row['key']:
            row.update(key=key,events=events,terminal=terminal)
    return sorted(groups.values(),key=lambda r:(r['key'],r['id']))


def merge_classes(ranked):
    """Union indexes only; never invert a compressed graph's generators."""
    groups = {}
    for classes in ranked:
        for row in classes:
            cid = row['id']
            entry = dict(direction=row['direction'],terminals=row['members'])
            if cid not in groups:
                groups[cid] = {**row,'origins':[entry]}
            else:
                old = groups[cid]
                origins = old['origins']+[entry]
                best = row if (row['key'],row['direction']) < (old['key'],old['direction']) else old
                groups[cid] = {**best,'origins':origins}
    return sorted(groups.values(),key=lambda r:(r['key'],r['id']))


def representative_metrics(classes, expected):
    rows = [dict(rank=i+1,score=row['key'][:3],reference_equivalent=row['id']==expected)
            for i,row in enumerate(classes)]
    return dict(class_count=len(rows),
        topk={str(k):any(r['reference_equivalent'] for r in rows[:k]) for k in (1,3,5,10)},
        **summarize_classes(rows))


def union_outcome(outcomes):
    if 'recovered' in outcomes:
        return 'recovered'
    if outcomes and all(o == 'not_recovered' for o in outcomes):
        return 'not_recovered'
    return 'unknown'




def summarize_classes(classes):
    """Event windows retain the best heavy/all-atom coverage tier."""
    best=classes[0]['score'][:3] if classes else None
    reference=next((c for c in classes if c['reference_equivalent']),None)
    windows={}
    for tolerance in (0,1,2,3,5):
        allowed=[c for c in classes if c['score'][:2]==best[:2]
                 and c['score'][2]<=best[2]+tolerance] if best else []
        windows[str(tolerance)]=dict(count=len(allowed),recovered=any(c['reference_equivalent'] for c in allowed))
    boundary=classes[min(5,len(classes))-1]['score'][:3] if classes else None
    tied=[c for c in classes if c['score'][:3]<=boundary] if boundary else []
    return dict(reference_class_rank=None if reference is None else reference['rank'],
        reference_event_gap=None if reference is None or reference['score'][:2]!=best[:2] else reference['score'][2]-best[2],
        top5_tie_inclusive=dict(count=len(tied),recovered=any(c['reference_equivalent'] for c in tied)),
        event_windows=windows)


def ranker(problem,heavy_only=False):
    endpoints=(problem.reactant,problem.product)
    edges=[[(a,b,float(e.wbo[a,b])) for a in range(e.atom_count) for b in range(a+1,e.atom_count)
        if e.wbo[a,b]>.2 and (not heavy_only or e.elements[a]!='H' and e.elements[b]!='H')]
        for e in endpoints]
    def score(mapping):
        mapping=dict(mapping)
        if heavy_only:mapping={r:p for r,p in mapping.items() if endpoints[0].elements[r]!='H' and endpoints[1].elements[p]!='H'}
        inverse={p:r for r,p in mapping.items()};broken=formed=changed=0
        for a,b,w in edges[0]:
            if a not in mapping or b not in mapping:broken+=int(a in mapping or b in mapping)
            elif endpoints[1].wbo[mapping[a],mapping[b]]<=.2:broken+=1
            elif abs(w-endpoints[1].wbo[mapping[a],mapping[b]])>.5:changed+=1
        for a,b,w in edges[1]:
            if a not in inverse or b not in inverse or endpoints[0].wbo[inverse.get(a,0),inverse.get(b,0)]<=.2:formed+=1
        heavy=sum(endpoints[1].elements[p]!='H' for p in mapping.values())
        return (-heavy,-len(mapping),broken+formed+changed,tuple(sorted(mapping.items()))),dict(broken=broken,formed=formed,order_changed=changed)
    return score
