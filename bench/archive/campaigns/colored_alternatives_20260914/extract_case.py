"""Read saved case 127 witnesses and supporting families; do not run AAM."""
import json,gzip,sys,hashlib,time,argparse
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--workspace',type=Path,required=True,help='Workspace with archived work/current-validation campaign')
parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2])
parser.add_argument('--package-src',type=Path,help='Optional already built package source tree')
args=parser.parse_args();ROOT=args.workspace.resolve();repo=args.repo.resolve()
if args.package_src:sys.path.insert(0,str(args.package_src.resolve()))
from mappa.artifacts import read_aam_checkpoint
from mappa.final_branches import FinalBranchCatalogue
from mappa.event_patterns import SignedEventIndex
from mappa.family_query import query_path
from rdkit import Chem
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
start=time.process_time();case=127
base=ROOT/f'work/current-validation/runs/coordinate/seed1/case{case}/R_to_P/cuts/aam.pkl.gz'
a=read_aam_checkpoint(base);problem=a.problem
cat=FinalBranchCatalogue(problem);cat.add_graph(a.graph,str(base))
# All three examples already have support in the base search archive.
source=repo/'reports/current_validation_20260912/coordinate-decoded.json.gz'
row=json.load(gzip.open(source))[case]
comparison_source=repo/'reports/current_validation_20260912/coordinate-slap-comparison.json.gz'
comparison=json.load(gzip.open(comparison_source))['rows'][case]
assert sha(base)==json.loads((base.parent.parent/'search.json').read_text())['archive_sha256']
idx=SignedEventIndex(problem)
patterns=[]
for p in row['patterns'].values():
 assert p['total']==4 and idx.describe(p['mapping'])['id']==p['id']
 fid=p['family'];f=cat.families[fid];m=dict(enumerate(p['mapping']))
 status,details=query_path(f.as_path(problem),problem,m,source_atoms=tuple(m),timeout_ms=5000)
 assert status=='recovered'
 groups=[key[1] for key,families in cat.branches.items() if fid in families]
 assert groups
 fragments=[list(pair[0]) for pair in groups[0]]
 patterns.append(dict(p,supporting_family=f.to_record(),source_fragments=fragments,query_status=status))
# A non-reacting methyl-H exchange within the first witness.
p=patterns[0];swapped=dict(enumerate(p['mapping']));swapped[9],swapped[10]=swapped[10],swapped[9]
status,detail=query_path(cat.families[p['family']].as_path(problem),problem,swapped,source_atoms=tuple(swapped),timeout_ms=5000)
assert status=='recovered' and idx.describe([swapped[i] for i in range(idx.n)])['id']==p['id']
def conventional(ep):
 m=Chem.RWMol()
 for i,z in enumerate(ep.elements):
  a=Chem.Atom(z);a.SetNoImplicit(True);a.SetAtomMapNum(i+1);m.AddAtom(a)
 for i in range(len(ep.elements)):
  for j in range(i+1,len(ep.elements)):
   w=ep.wbo[i,j]
   if w>=.5:m.AddBond(i,j,Chem.BondType.DOUBLE if w>=1.4 else Chem.BondType.SINGLE)
 m=m.GetMol();Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
 assert m.GetNumAtoms()==23 and all(a.GetNumRadicalElectrons()==0 for a in m.GetAtoms())
 return Chem.MolToSmiles(m)
result=dict(case=case,collection='140 coordinate reactions; mapping-unverified development collection',
 source_sha256={str(base.relative_to(ROOT)):sha(base),str(source.relative_to(repo)):sha(source),str(comparison_source.relative_to(repo)):sha(comparison_source)},
 input=json.loads((ROOT/f'work/full140_inputs/{case}/input.json').read_text()),
 reactant_smiles=conventional(problem.reactant),product_smiles=conventional(problem.product),
 depiction='Bond-order threshold 1.4 for conventional single/double depictions, adjacency 0.5. Colors/counts derive independently from original WBO matrices at event thresholds 0.5/0.3.',
 patterns=patterns,minimum=4,window=row['window'],saved_window_complete=row['complete_saved_window'],
 comparator_minimum_ids={k:list(v['patterns']) for k,v in comparison['methods'].items()},
 shuffle=dict(source_atoms=[9,10],mapping=swapped,status='allowed',same_event_id=p['id']),
 search_runs=0,decode_runs=0,cpu_seconds=time.process_time()-start)
out=Path(__file__).resolve().parent/'coordinate_case127.json';out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'case':case,'seconds':result['cpu_seconds'],'patterns':[(p['total'],p['source_fragments']) for p in patterns],'R':result['reactant_smiles'],'P':result['product_smiles'],'shuffle':status},indent=2))
