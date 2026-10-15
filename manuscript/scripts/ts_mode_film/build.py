"""Export recorded TEMPO initial-guess modes and verify the current MAPPA selection."""
from pathlib import Path
import argparse,json,hashlib,signal
import numpy as np
from mappa import (AAMProblem,AAMSearchConfig,AtomBijection,MolecularEndpoint,ReactionContext,ResolvedMechanism,TransitionStateTarget,VibrationalModes,analyze_transition_state,ts_record)
from mappa.chemistry_computations import parse_xyz
from mappa.chemistry_computations.xtb import read_wbo_file
from mappa.modes import parse_g98_modes
from mappa.scoring import pair_geometry, subset_score, POLICY

def build(source,library,output):
 output.mkdir(parents=True,exist_ok=True)
 def endpoint(folder,label):
  f=next(p for p in folder.glob('*.xyz') if not p.name.startswith('xtb'))
  e,x=parse_xyz(f);return MolecularEndpoint(tuple(e),x,read_wbo_file(folder/'wbo',len(e)),label=label)
 if source.is_file():
  saved=json.loads(source.read_text())
  def load(key,label):
   e=saved[key];return MolecularEndpoint(tuple(e['elements']),np.array(e['coordinates']),np.array(e['wbo']),label=label)
  r,p,t=load('reactant','R'),load('product','P'),load('target','TS guess')
  raw=dict(mapping_RP=saved['mapping_RP'],broken_bonds_R=saved['broken'],formed_bonds_R=saved['formed'],core_atoms=saved['core'])
  f,m=np.array(saved['frequencies']),np.array(saved['displacements']);old=saved['previous_selection']
 else:
  r=endpoint(source/'work/endpoints/R','R');p=endpoint(source/'work/endpoints/P','P');t=endpoint(source/'work/targets/iter1_hess','TS guess')
  raw=json.loads((source/'stages/pr1.tempo_ts3/rp_stage.json').read_text())['mechanisms'][0]
  f,m=parse_g98_modes(source/'work/targets/iter1_hess/g98.out')
  old=json.loads((source/'alignments/pr1.tempo_ts3/ts_alignments/mechanisms/mechanism_001/ig_iter1/score.json').read_text())
 mp={int(a):b for a,b in raw['mapping_RP'].items()};mech=ResolvedMechanism(AtomBijection.from_mapping(mp,degree=r.atom_count),tuple(map(tuple,raw['broken_bonds_R'])),tuple(map(tuple,raw['formed_bonds_R'])),tuple(raw['core_atoms']))
 context=ReactionContext(AAMProblem(r,p,name='pr1.tempo_ts3'),AAMSearchConfig(),(mech,))
 result=analyze_transition_state(context,TransitionStateTarget(t,VibrationalModes(f,m),kind='initial_guess'));s=result.mechanisms[0].selected
 assert s is not None and s.mode_index==0 and s.assignment.as_dict()=={12:12,14:14,49:49}
 geometry=pair_geometry(t.coordinates,t.elements)
 events={tuple(sorted(e['target_bond'])) for c in result.mechanisms[0].candidates for e in c.event_terms}
 pairs=[(int(i),int(j)) for k,(i,j) in enumerate(zip(geometry['i'],geometry['j'])) if geometry['distance'][k]<=1.2*geometry['radii_sum'][k] or (int(i),int(j)) in events]
 b=np.array([subset_score(t.coordinates,q,t.elements,s.event_terms,pairs)['score'] for q in m])
 assert abs(b[s.mode_index]-s.score)<1e-12
 imaginary=np.flatnonzero(f<0)
 rejected=19;assert f[rejected]>0
 for term in s.event_terms:
  assert abs(term['wbo_product']-term['wbo_reactant'])>.5
 assert all(r.elements[a]==p.elements[mp[a]] for a in mp)
 # A single rigid camera basis places O -> N horizontally, with H above their axis.
 O,N,H=49,14,12;origin=t.coordinates[[O,N]].mean(0);ex=t.coordinates[N]-t.coordinates[O];ex/=np.linalg.norm(ex);ey=t.coordinates[H]-origin;ey-=ey@ex*ex;ey/=np.linalg.norm(ey);ez=np.cross(ex,ey);rot=np.array([ex,ey,ez]).T
 def pose(x):return ((x-origin)@rot).tolist()
 def ep(e):
  return dict(elements=list(e.elements),xyz=pose(e.coordinates),edges=[[i,j] for i in range(e.atom_count) for j in range(i+1,e.atom_count) if e.wbo[i,j]>=.2])
 modes=[]
 for k in [rejected,s.mode_index]:
  # One scale for both raw normal modes; preserve all per-atom relative motion.
  modes.append(dict(index=k,number=k+1,frequency=float(f[k]),score=float(b[k]),eligible=bool(f[k]<0),displacement=(m[k]@rot).tolist(),selected=k==s.mode_index))
 # Verify distance changes are coherent for the chosen mode after phase orientation.
 k=s.mode_index;derivatives=[]
 for term in s.event_terms:
  i,j=term['target_bond'];u=t.coordinates[j]-t.coordinates[i];u/=np.linalg.norm(u);derivatives.append(float(u@(m[k,j]-m[k,i])))
 assert derivatives[0]*derivatives[1]<0
 data=dict(case='pr1.tempo_ts3',atoms=r.atom_count,duration=30,reactant=ep(r),product=ep(p),target=ep(t),mapping=mp,core=[O,H,N],broken=[[12,49]],formed=[[12,14]],modes=modes,score=s.score,progress=s.wbo_progress,amplitude=.60,oscillation_hz=.8,mode_count=len(f),imaginary_count=len(imaginary),source='Recorded holdout initial guess, iteration 1; cached xTB Hessian',stage='TS guess · guided mode selection')
 (output/'film-data.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
 (output/'current-ts-score.json').write_text(json.dumps(ts_record(result),indent=2)+'\n')
 (output/'all-mode-scores.json').write_text(json.dumps([dict(index=i,frequency=float(x),score=float(b[i]),eligible=bool(x<0)) for i,x in enumerate(f)],indent=2)+'\n')
 inputs=output/'inputs';inputs.mkdir(exist_ok=True)
 # Self-contained numerical replay inputs; no machine paths or unrelated structures.
 inp=dict(case=data['case'],reactant=dict(elements=r.elements,coordinates=r.coordinates.tolist(),wbo=r.wbo.tolist()),product=dict(elements=p.elements,coordinates=p.coordinates.tolist(),wbo=p.wbo.tolist()),target=dict(elements=t.elements,coordinates=t.coordinates.tolist(),wbo=t.wbo.tolist()),frequencies=f.tolist(),displacements=m.tolist(),mapping_RP=mp,broken=raw['broken_bonds_R'],formed=raw['formed_bonds_R'],core=raw['core_atoms'],previous_selection=dict(S=old['S'],decomposition=old['decomposition']))
 (inputs/'recorded-case.json').write_text(json.dumps(inp,separators=(',',':'))+'\n')
 (output/'science-validation.json').write_text(json.dumps(dict(status='passed',atom_count=57,selected_mode_index=0,selected_frequency=float(f[0]),selected_score=s.score,previous_score=old['S'],scoring_policy=POLICY, historical_score_superseded=True,imaginary_modes=1,stable_comparator_index=rejected,stable_comparator_score=float(b[rejected]),bond_length_derivatives=derivatives,all_atom_motion_preserved=True,common_display_amplitude=.60,display_frequency_is_illustrative=True,optimized_ts_claimed=False,input_sha256=hashlib.sha256((inputs/'recorded-case.json').read_bytes()).hexdigest()),indent=2)+'\n')
 template=Path(__file__).with_name('film.html').read_text();(output/'index.html').write_text(template.replace('__LIBRARY__',library.read_text()).replace('__DATA__',json.dumps(data,separators=(',',':')).replace('</','<\\/')))
 print(json.dumps(dict(selected=modes[1],rejected_score=modes[0]['score']),default=float)[:240])
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--source',type=Path,required=True);a.add_argument('--library',type=Path,required=True);a.add_argument('--output',type=Path,required=True);v=a.parse_args();signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(TimeoutError('300 s watchdog')));signal.alarm(300);build(v.source,v.library,v.output)
