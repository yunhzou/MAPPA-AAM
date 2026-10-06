from pathlib import Path
import sys,json,os
PACKAGE=Path(__file__).resolve().parent;S=Path(os.environ.get('MAPPA_EXPERIMENT_WORK',PACKAGE/'work')).resolve();sys.path.insert(0,str(S/'engine/src'))
import numpy as np
from mappa.frag import build_graph
from mappa.fragment import match_fragment,FragmentMatchContext,FragmentMatchConfig
from mappa.fragment_choices import FragmentChoices

def test_mapped_prefixes():
 for n in [5,8]:
  w=np.zeros((n,n))
  for i in range(n-1):w[i,i+1]=w[i+1,i]=1
  # Distinct element labels make the expected growth order unambiguous.
  g=build_graph([f'E{i}' for i in range(n)],w,.2)
  initial={0:0,1:1};ctx=FragmentMatchContext(initial,{0:0,1:0})
  cfg=FragmentMatchConfig(iso_tolerance=1,branch_limit=2000,allow_mapped_seed=True)
  expected=match_fragment(g,g,seed=0,context=ctx,config=cfg)
  session=FragmentChoices(g,g,seed=0,context=ctx,config=cfg)
  assert session.normal.matches==expected.matches and not session.normal.capped
  seen=[]
  for c in session.closures:
   for match in session.close(c).matches:
    merged=dict(initial);merged.update(dict(match));assert all(merged[r]==p for r,p in initial.items())
    assert set(merged)==set(range(max(merged)+1));seen.append(len(merged))
  assert set(range(2,n))<=set(seen)
 print('Mapped-seed observer preserves the greedy result and every chain prefix.')

if __name__=='__main__':test_mapped_prefixes()
