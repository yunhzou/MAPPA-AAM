"""Golden reference and selected current MAPPA alternatives; no search at build time."""
from pathlib import Path
import json,os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Circle
from rdkit import Chem
from rdkit.Chem import rdDepictor
from build_molecule_figure import molecule,bond_event_marker,INK,MUTED,TEAL,RED,BLUE
R='[C:1]([O:2][H:6])([H:3])([H:4])[H:5]'
A='[C:1](=[O:2])([H:3])[H:4]'
B='[C:1](=[O:2])([H:3])[H:6]'

def check_example(man):
 d=json.loads((man/'evidence/multicandidate_example.json').read_text())
 assert d['complete_saved_families'] and [c['total'] for c in d['candidates']]==[4,6]
 assert d['queries']=={'same_class_H3_H4_swap':'allowed','other_class_preserve_events':'forbidden','other_class_any_saved_events':'allowed'}
 wr=np.array(d['wbo_R']);wp=np.array(d['wbo_P'])
 expected=[[(0,1,1),(0,4,-1),(1,5,-1),(4,5,1)],[(0,1,1),(0,3,-1),(0,4,-1),(0,5,1),(1,5,-1),(3,4,1)]]
 from rdkit import Chem
 def bonds(s):
  p=Chem.SmilesParserParams();p.removeHs=False;m=Chem.MolFromSmiles(s,p);assert m
  return {tuple(sorted((b.GetBeginAtom().GetAtomMapNum()-1,b.GetEndAtom().GetAtomMapNum()-1))):b.GetBondTypeAsDouble() for b in m.GetBonds()}
 r=bonds(R)
 assert len(r)==5 and all(wr[i,j]==v for (i,j),v in r.items())
 for m,expected_events,smiles in zip(d['figure_witnesses_zero_based'],expected,[A+'.[H:5][H:6]',B+'.[H:4][H:5]']):
  m={int(k):v for k,v in m.items()}
  delta=wp[np.ix_([m[i] for i in range(6)],[m[i] for i in range(6)])]-wr
  events=[(i,j,int(np.sign(delta[i,j]))) for i in range(6) for j in range(i+1,6) if abs(delta[i,j])>=.5]
  assert events==expected_events
  b=bonds(smiles)
  assert len(b)==4
  assert all(b.get((i,j),0)==wp[m[i],m[j]] for i in range(6) for j in range(i+1,6))
 counts=json.loads((man/'evidence/output_multiplicity.json').read_text())
 assert all(m['cases']==1851 for m in counts['methods'])
 assert all(m['multiple_record_cases']==0 for m in counts['methods'] if not m['method'].startswith('slap'))
 assert all(m['multiple_record_cases']>0 for m in counts['methods'] if m['method'].startswith('slap'))
 return d


def check_case(man):
 d=json.loads((man/'evidence/golden_case9.json').read_text())
 proof=json.loads((man/'evidence/golden_case9_verification.json').read_text())
 assert d['case']==9 and proof['original_rdf_annotation_verified']
 assert d['config']['seed_count']==1 and d['config']['branch_limit']==2000 and d['config']['sweep_cuts']
 assert d['patterns'][0]['exact_reference_heavy']
 ref={int(k):v for k,v in d['reference_mapping'].items()}
 first={int(k):v for k,v in d['patterns'][0]['mapping'].items()}
 assert all(first[k]==v for k,v in ref.items())
 assert [p['total'] for p in d['patterns']]==[9,8,9]
 assert [p['heavy_changes'] for p in proof['rows']]==[5,6,7]
 assert all(p['attains_H_lower_bound'] for p in proof['rows'])
 from rdkit import Chem
 wr=np.array(d['input']['reactant']['wbo']);wp=np.array(d['input']['product']['wbo'])
 for p in d['patterns']:
  m={int(k):v for k,v in p['mapping'].items()};v=[m[i] for i in range(33)]
  assert sorted(v)==list(range(33))
  assert all(d['input']['reactant']['elements'][i]==d['input']['product']['elements'][j] for i,j in m.items())
  assert sorted(a for f in p['source_fragments'] for a in f)==list(range(33))
  delta=wp[np.ix_(v,v)]-wr
  hydrogens_r=[i for i,e in enumerate(d['input']['reactant']['elements']) if e=='H']
  hydrogens_p=[i for i,e in enumerate(d['input']['product']['elements']) if e=='H']
  heavy=[i for i in range(33) if i not in hydrogens_r]
  for matrix,hs in [(wr,hydrogens_r),(wp,hydrogens_p)]:
   assert all(sum(matrix[i,j]>.5 for j in range(33))==1 for i in hs)
   assert all(matrix[i,j]==0 for i in hs for j in hs)
  lower=sum(abs(sum(wr[i,h]>.5 for h in hydrogens_r)-sum(wp[m[i],h]>.5 for h in hydrogens_p)) for i in heavy)
  assert lower==p['total']-sum(map(len,p['heavy_events'].values()))
  assert p['events']=={key:[[i,j] for i in range(33) for j in range(i+1,33) if sign*delta[i,j]>=.5] for key,sign in [('broken',-1),('formed',1)]}
 for side,w in [('reactant',wr),('product',wp)]:
  params=Chem.SmilesParserParams();params.removeHs=False;mol=Chem.MolFromSmiles(d[side+'_smiles'],params)
  matrix=np.zeros((33,33))
  for b in mol.GetBonds():
   i=b.GetBeginAtom().GetAtomMapNum()-1;j=b.GetEndAtom().GetAtomMapNum()-1
   matrix[i,j]=matrix[j,i]=b.GetBondTypeAsDouble()
  assert np.array_equal(matrix,w)
 return d

def build(man):
 D=check_case(man);labels=D['original_labels'][0]
 PALETTE=[(.70,.87,.81),(.98,.81,.49),(.81,.68,.92),(.66,.81,.95),(.96,.73,.78),(.81,.83,.66)]
 params=Chem.SmilesParserParams();params.removeHs=False
 base=[Chem.RemoveHs(Chem.MolFromSmiles(D[side+'_smiles'],params)) for side in ['reactant','product']]
 coords=[]
 for mol in base:
  rdDepictor.Compute2DCoords(mol)
  c={a.GetAtomMapNum():tuple(mol.GetConformer().GetAtomPosition(a.GetIdx()))[:2] for a in mol.GetAtoms()}
  for frag in Chem.GetMolFrags(mol):
   keys=[mol.GetAtomWithIdx(i).GetAtomMapNum() for i in frag]
   cx=(min(c[k][0] for k in keys)+max(c[k][0] for k in keys))/2
   cy=(min(c[k][1] for k in keys)+max(c[k][1] for k in keys))/2
   target_x=(8 if len(frag)==6 else 0) if len(coords)==0 else (6.5 if len(frag)==1 else 0)
   for key in keys:c[key]=(c[key][0]-cx+target_x,c[key][1]-cy)
  coords.append(c)
 fig=plt.figure(figsize=(11.8,9.0));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1416),ylim=(1080,0));ax.axis('off')
 def t(x,y,s,size=11,color=INK,weight='normal',ha='left'):
  ax.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va='center',linespacing=1.3,zorder=9)
 def box(x,y,w,h,color='#F5F7FB',edge=None):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=12',facecolor=color,edgecolor=edge or color,lw=.8,zorder=0))
 def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,color=INK,lw=1.8,zorder=6))
 box(12,12,1392,1056)
 t(32,41,'Golden case 9: the reference is one of several recovered mappings',15,weight='bold')
 t(32,73,'MAPPA search · 1 seed (default) + cut sweep · cap 2,000 · iso tolerance 1 · event threshold 0.5',10.5,MUTED)
 for x,c in zip([44,69,94],PALETTE):ax.add_patch(Circle((x,112),10,facecolor=c,edgecolor='none'))
 t(116,112,'Saved fragment groups',10,MUTED)
 for x,sign,txt in [(526,-1,'Break / weaken'),(972,1,'Form / strengthen')]:
  ax.plot([x,x+64],[112,112],color='black',lw=1.1)
  bond_event_marker(ax,(x,112),(x+64,112),sign,size=9)
  t(x+79,112,txt,10,INK)
 row_titles=['A   Ground truth (Golden reference) — recovered by MAPPA','B   MAPPA alternative — different oxygen retained','C   MAPPA alternative — oxygen and carbon correspondence change']
 wr=np.array(D['input']['reactant']['wbo']);wp=np.array(D['input']['product']['wbo'])
 for k,p in enumerate(D['patterns']):
  y=142+k*275;box(28,y,1359,263,'white','#C7DFD5' if k==0 else '#DCE3E9')
  t(47,y+24,row_titles[k],12,TEAL if k==0 else INK,'bold')
  heavy=sum(map(len,p['heavy_events'].values()));hyd=p['total']-heavy
  t(1358,y+24,f"{heavy} heavy-atom changes",10.5,TEAL,'bold','right')
  m={int(i):j for i,j in p['mapping'].items()};inverse={j:i for i,j in m.items()}
  assert sorted(m)==list(range(33)) and sorted(m.values())==list(range(33))
  v=[m[i] for i in range(33)];delta=wp[np.ix_(v,v)]-wr
  assert p['events']=={key:[[i,j] for i in range(33) for j in range(i+1,33) if sign*delta[i,j]>=.5] for key,sign in [('broken',-1),('formed',1)]}
  owners={labels[i]:PALETTE[q%len(PALETTE)] for q,f in enumerate(p['source_fragments']) for i in f if i<15}
  for side,(x,w) in enumerate([(48,570),(765,590)]):
   mol=Chem.Mol(base[side]);co={}
   for at in mol.GetAtoms():
    original=at.GetAtomMapNum();src=original-1 if side==0 else inverse[original-1]
    label=labels[src];at.SetAtomMapNum(label);co[label]=coords[side][original]
   ev=p['heavy_events']['broken' if side==0 else 'formed'];changes={tuple(labels[i] for i in pair):(-1 if side==0 else 1) for pair in ev}
   positions=molecule(ax,Chem.MolToSmiles(mol),x,y+48,w,183,owners=owners,font=26,align_acetyl=False,bond_changes=changes,coordinates=co)
   if side==0:
    px=(max(positions[i][0] for i in range(1,10))+min(positions[i][0] for i in range(10,16)))/2
   else:
    water_label=labels[inverse[14]]
    px=(max(xy[0] for label,xy in positions.items() if label!=water_label)+positions[water_label][0])/2
   t(px-(14 if side else 0),y+139,'+',15,MUTED,ha='center')
  arrow((660,y+139),(727,y+139))
  t(47,y+242, ['Reference: alcohol O₁₅ enters the ring; ketone O₃ leaves in water.','Alternative: ketone O₃ enters the ring; alcohol O₁₅ leaves in water.','Alternative: O₃ enters the ring, with a different assignment of C₁₃ and C₁₄.'][k],10,MUTED)
  t(1358,y+242,f"{p['total']} total events = {heavy} heavy + {hyd} involving H",9.5,MUTED,ha='right')
 t(40,1000,'Retaining alternatives matters: a lower total event count (B: 8) need not select the annotated mapping (A: 9).',11,TEAL,'bold')
 t(40,1031,'Ground truth refers to Golden’s heavy-atom annotation. H identities are unannotated; totals use each MAPPA witness.',9.7,MUTED)
 t(40,1052,'Numbers are original Golden labels. H atoms are implicit. Selected witnesses do not imply exhaustive or minimum-event decoding.',9.7,MUTED)
 for ext in ['pdf','svg','png']:fig.savefig(man/'figs'/f'fig4_alternatives.{ext}',facecolor='white',bbox_inches='tight',pad_inches=.025,dpi=220)
 plt.close(fig)
 print('Figure generated')


if __name__=="__main__":build(Path(__file__).resolve().parents[1])
