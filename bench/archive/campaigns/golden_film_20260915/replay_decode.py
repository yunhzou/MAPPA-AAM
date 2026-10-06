"""Verify the Golden film's saved competition and complete event-window decode.

Run from reports/golden_film_20260915 in an installed repository checkout.
No new mapping search or bijection enumeration is performed.
"""
from pathlib import Path
import gzip,json,sys,time
S=Path(__file__).resolve().parent
sys.path.insert(0,str(S.parents[1]/'src'))
from mappa.artifacts import aam_from_record
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex,extract_path_events
from mappa.competition import takeover_plan

read=lambda name:json.loads(gzip.decompress((S/name).read_bytes()))
a=aam_from_record(read('baseline.json.gz'))
r=aam_from_record(read('takeover3.json.gz'))
d=json.loads((S/'film-source.json').read_text());proof=d['competition']['proof']
assert proof['parent']==29 and proof['repair_id']=='takeover3'
old=dict(a.graph.states[29].mapping);owners=dict(a.graph.states[29].islands)
plan=takeover_plan(old,owners,owners[proof['seed']],dict(proof['B_placement']),allow_b_relocation=True)
for atom in proof['released_dependents']:
    plan['anchors'].pop(atom,None)
assert plan['anchors']==dict(r.graph.contexts[0].anchors)
assert sorted(set(old)-plan['anchors'].keys())==proof['holes']==[17,19]
assert dict(r.graph.states[3].mapping)=={int(k):v for k,v in d['paths'][2]['mapping'].items()}
cat=FinalBranchCatalogue.from_record(a.problem,read('catalogue.json.gz'))
index=SignedEventIndex(a.problem,threshold=.5,metal_threshold=.3)
start=time.process_time();ids=set()
for i,family in enumerate(cat.families):
    family.validate_representative(a.problem)
    result=extract_path_events(family.as_path(a.problem),a.problem,index,max_events=6,seconds=5,max_patterns=None)
    assert result['complete'],(i,result['reason'])
    for p in result['patterns']:
        assert index.describe(p['mapping'])['id']==p['id']
        ids.add(p['id'])
expected=json.loads((S/'decode.json').read_text())
assert ids==set(expected['patterns']) and len(ids)==9
assert cat.branch_count==58 and len(cat.families)==118
for candidate in d['decoded_candidates']:
    family=cat.families[candidate['family']]
    result=extract_path_events(family.as_path(a.problem),a.problem,index,max_events=6,seconds=5,max_patterns=None)
    vector=[candidate['mapping'][str(i)] for i in range(index.n)]
    assert any(p['mapping']==vector and p['id']==candidate['event_class'] for p in result['patterns'])
assert d['paths'][1]['event_class']==d['paths'][2]['event_class']
assert d['paths'][0]['event_class']!=d['paths'][1]['event_class']
# Use the same charged, heavy-atom graph criterion as the Golden evaluation.
sys.path.insert(0,str(S.parents[1]/'bench'))
from golden_evaluation import colored_graph,project
import pynauty
reference=json.loads((S/'reference.json').read_text())
features=reference['features']
certificate=lambda m:pynauty.certificate(colored_graph(features,project(m,features)))
truth=certificate(dict(reference['mapping']))
for c in d['decoded_candidates']:
    assert (certificate({int(k):v for k,v in c['mapping'].items()})==truth)==c['reference_equivalent']
print(json.dumps(dict(status='passed',branches=cat.branch_count,families=len(cat.families),event_classes=len(ids),shown=2,competition_released=proof['holes'],cpu_seconds=time.process_time()-start)))
