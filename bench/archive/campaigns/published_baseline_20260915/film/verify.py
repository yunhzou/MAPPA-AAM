"""Verify baseline-only film witnesses and reference labels; no search."""
from pathlib import Path
import sys,json,gzip
S=Path(__file__).resolve().parent;R=S.parents[2];sys.path.insert(0,str(R/'src'))
from mappa.artifacts import aam_from_record
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex
from mappa.family_query import query_path
read=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
a=aam_from_record(read(S/'baseline.json.gz'));cat=FinalBranchCatalogue.from_record(a.problem,read(S/'catalogue.json.gz'));d=read(S/'film-source.json');idx=SignedEventIndex(a.problem);assert 'competition' not in d
assert len(d['paths'])==2 and {p['terminal'] for p in d['paths']}=={29,30}
assert len({p['event_class'] for p in d['paths']})==2
assert all(p['event_class']==d['decoded_candidates'][p['class_index']]['event_class'] for p in d['paths'])
for p in d['paths']:
 m={int(k):v for k,v in p['mapping'].items()};assert m==dict(a.graph.states[p['terminal']].mapping)
 assert idx.describe([m[i] for i in range(idx.n)])['id']==p['event_class']
for c in d['decoded_candidates']:
 m={int(k):v for k,v in c['mapping'].items()};assert query_path(cat.families[c['family']].as_path(a.problem),a.problem,m,source_atoms=tuple(range(idx.n)))[0]=='recovered'
 raw=idx.describe([m[i] for i in range(idx.n)]);assert raw['id']==c['event_class']
 assert {(e['kind'],tuple(e['r'])) for e in c['events']}=={(k,tuple(p)) for k,ps in raw['events'].items() for p in ps}
sys.path.insert(0,str(R/'bench'));from golden_evaluation import colored_graph,project
import pynauty
ref=read(S/'reference.json');cert=lambda m:pynauty.certificate(colored_graph(ref['features'],project(m,ref['features'])))
for c in d['decoded_candidates']:assert (cert({int(k):v for k,v in c['mapping'].items()})==cert(dict(ref['mapping'])))==c['reference_equivalent']
assert cat.branch_count==44 and len(cat.families)==96
assert len(read(S/'decode.json')['patterns'])==9 and len({c['event_class'] for c in d['decoded_candidates']})==2
print('Passed: two baseline witnesses, two distinct displayed classes, reference label, 44 branches / 96 families.')
