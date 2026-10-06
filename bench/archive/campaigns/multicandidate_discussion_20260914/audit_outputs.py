"""Small current-API reproduction plus read-only archived-output audit."""
import sys,json,gzip,hashlib,time,argparse
from pathlib import Path
from collections import Counter
from dataclasses import asdict
parser=argparse.ArgumentParser()
parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2])
parser.add_argument('--package-src',type=Path,help='Optional built source tree containing mappa')
parser.add_argument('--archive-dir',type=Path)
args=parser.parse_args();repo=args.repo.resolve()
archive=args.archive_dir or repo/'reports/golden_competitor_recheck_20260913'
sys.path.insert(0,str(repo/'docs'))
if args.package_src:
    package=args.package_src.resolve()
    for f in (repo/'src/mappa').rglob('*.py'):
        assert f.read_bytes()==(package/'mappa'/f.relative_to(repo/'src/mappa')).read_bytes(),f
    sys.path.insert(0,str(package))
from mappa import AAMProblem,AAMSearchConfig,search_aam
from mappa.postprocessing import EventDecodeConfig,decode_events
from notebook_helpers import endpoint
out=Path(__file__).resolve().parent
p=AAMProblem(endpoint('CO','Methanol'),endpoint('C=O.[H][H]','Formaldehyde + H2'))
config=AAMSearchConfig(iso_tolerance=1,seed_count=1,branch_limit=2000,random_seed=42)
t=time.process_time();w=time.perf_counter()
a=search_aam(p,config,workers=1)
d=decode_events(a,EventDecodeConfig(threshold=.5,metal_threshold=.3))
assert d.complete and [c.total for c in d.candidates]==[4,6]
# Fix source-index representatives for the figure, irrespective of solver witness order.
witnesses=[{0:0,1:1,2:4,3:5,4:3,5:2},{0:0,1:1,2:5,3:3,4:2,5:4}]
checks=[]
for c,m in zip(d.candidates,witnesses):
 q=d.query(c,m)
 assert q.status=='allowed'
 checks.append(q.status)
c=d.candidates[0];m=witnesses[0];swapped={**m,2:m[3],3:m[2]}
queries={'same_class_H3_H4_swap':d.query(c,swapped).status,
 'other_class_preserve_events':d.query(c,witnesses[1]).status,
 'other_class_any_saved_events':d.query(c,witnesses[1],same_events=False).status}
assert queries=={'same_class_H3_H4_swap':'allowed','other_class_preserve_events':'forbidden','other_class_any_saved_events':'allowed'}
proof={'description':'Six-atom notebook reaction, current library, one direction; formal bond orders; not a pathway validation',
 'config':asdict(config),'decode_config':asdict(d.config),'complete_saved_families':d.complete,
 'family_count':len(d.catalogue.families),'candidates':[asdict(c) for c in d.candidates],
 'figure_witnesses_zero_based':witnesses,'witness_membership':checks,'queries':queries,
 'elements_R':list(p.reactant.elements),'elements_P':list(p.product.elements),
 'wbo_R':p.reactant.wbo.tolist(),'wbo_P':p.product.wbo.tolist(),
 'cpu_seconds':time.process_time()-t,'wall_seconds':time.perf_counter()-w,
 'source_sha256':{str(f.relative_to(repo)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [repo/'src/mappa/aam.py',repo/'src/mappa/postprocessing.py',repo/'docs/notebook_helpers.py']}}
(out/'multicandidate_example.json').write_text(json.dumps(proof,indent=2)+'\n')
methods=[]
for name in ['rxnmapper','localmapper','indigo','chython','rdt','slap_binary','slap_weighted']:
 f=archive/f'{name}-rescored.json.gz'
 records=json.load(gzip.open(f))['records']
 counts=Counter(r.get('candidates',0) for r in records)
 methods.append({'method':name,'cases':len(records),'candidate_count_histogram':dict(sorted(counts.items())),
   'multiple_record_cases':sum(v for k,v in counts.items() if k>1),'archive_sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
 if not name.startswith('slap'):assert max(counts)==1
summary={'scope':'Archived evaluated interfaces, 1851 reactions each; counts are returned explicit candidate records, not distinct event classes or expanded symmetry mappings. No comparator rerun.',
 'methods':methods,'harness':'bench/golden_competitors.py','harness_sha256':hashlib.sha256((repo/'bench/golden_competitors.py').read_bytes()).hexdigest()}
(out/'output_multiplicity.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'example_seconds':proof['cpu_seconds'],'family_count':proof['family_count'],'queries':queries,'method_counts':methods},indent=2))
