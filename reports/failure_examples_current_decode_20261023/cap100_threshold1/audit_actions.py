"""Verify saved constructive symmetry certificates without solving or enumerating mappings."""
import argparse,gzip,json,time
from pathlib import Path
from decode_shard import load
from mappa.family_query import group_factors


def program(family,problem):
    path=family.as_path(problem);steps=[];n=problem.atom_count
    for eid in reversed(path.transitions):
        placement=path.graph.fragment_placement(eid)
        if placement is None:continue
        for domain in placement.symmetry_domains:
            if domain.source=='exact_automorph_group' or len(domain.p_atoms)<2:continue
            pool=set(domain.p_atoms)
            locked=dict(path.graph.states[path.graph.transitions[eid].source].mapping)
            block_atoms={r for block in placement.symmetry_domains for r in block.r_atoms if block.source!='exact_automorph_group'}
            witness=dict(placement.representative_assignments)
            fixed=set(locked.values())|{witness[r] for r in placement.exact_fixed}|{p for r,p in witness.items() if r not in block_atoms}
            steps.append(('domain',eid,pool,fixed&pool))
        generators=tuple(g.images for g in placement.target_generators)
        for factor in group_factors(generators):
            if any(g[p]!=p for p in range(min(n,len(factor[0]))) for g in factor):
                steps.append(('group',eid,set(factor),set()))
    return steps


def verify(family,problem,pattern,steps=None):
    actual=dict(family.mapping);mapping=dict(enumerate(pattern['mapping']))
    if not pattern['actions']:
        assert actual==mapping
    else:
        steps=program(family,problem) if steps is None else steps
        assert len(pattern['actions'])==len(steps)
        for action,(kind,eid,allowed,fixed) in zip(pattern['actions'],steps):
            assert action['kind']==kind and action['transition']==eid
            pairs=action['permutation'];permutation=dict(pairs)
            assert len(permutation)==len(pairs)
            if kind=='domain':
                assert set(permutation)==allowed and set(permutation.values())==allowed
                assert all(permutation[p]==p for p in fixed)
            else:
                degree=len(next(iter(allowed)))
                assert set(permutation)==set(range(degree))
                assert tuple(permutation[p] for p in range(degree)) in allowed
            actual={r:permutation.get(p,p) for r,p in actual.items()}
        assert actual==mapping
    assert sorted(mapping)==list(range(problem.atom_count))
    assert sorted(mapping.values())==list(range(problem.atom_count))
    assert all(problem.reactant.elements[r]==problem.product.elements[p] for r,p in mapping.items())
    floor,iso=family.policy
    assert all(problem.product.wbo[mapping[a],mapping[b]]>=floor and abs(problem.reactant.wbo[a,b]-problem.product.wbo[mapping[a],mapping[b]])<=iso+1e-9 for a,b in family.required_edges)


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);a=p.parse_args();start=time.perf_counter()
    out=a.run/'results/search_memory_partitioned';cat=load(out,a.inputs);d=json.load(gzip.open(out/'decoded.json.gz','rt'));wanted={};identity=0;checked=set();actions=0
    for c in d['candidates']:
        fid=c['family_ids'][0];mapping={int(k):v for k,v in c['mapping'].items()}
        if mapping==dict(cat.families[fid].mapping):
            verify(cat.families[fid],cat.problem,dict(mapping=[mapping[i] for i in range(cat.problem.atom_count)],actions=[]));checked.add(c['id']);identity+=1
        else:wanted.setdefault(fid,{})[c['id']]=mapping
    for path in sorted((a.run/'shards').glob('shard-*.jsonl')):
        with path.open() as f:
            for line in f:
                row=json.loads(line);fid=row['family_id']
                if fid not in wanted:continue
                family=cat.families[fid];steps=program(family,cat.problem)
                for pattern in row['patterns']:
                    cid=pattern['id']
                    if cid not in wanted[fid]:continue
                    assert dict(enumerate(pattern['mapping']))==wanted[fid][cid]
                    assert cid not in checked
                    verify(family,cat.problem,pattern,steps)
                    checked.add(cid);actions+=1
    assert checked=={c['id'] for c in d['candidates']}
    result=dict(status='passed',witnesses=len(checked),identity_action_witnesses=identity,constructive_action_certificates=actions,minimum_witnesses_verified=sum(c['total']==d['observed_minimum'] for c in d['candidates']),checks=['permutation domains and fixed atoms','exact membership in recorded group factors','action order and transition identity','exact replay to emitted full mapping','element preservation and injectivity','every required fragment edge and WBO tolerance'],seconds=time.perf_counter()-start)
    (out/'constructive-membership-audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
