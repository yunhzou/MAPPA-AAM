from pathlib import Path
import sys,os,importlib.util,json
PACKAGE=Path(__file__).resolve().parent;S=Path(os.environ.get('MAPPA_EXPERIMENT_WORK', PACKAGE/'work')).resolve();sys.path.insert(0,str(S/'engine/src'))
from mappa.artifacts import read_aam_checkpoint
from mappa.competition import compete_fragments,CompetitionConfig
name='mappa._control_competition';spec=importlib.util.spec_from_file_location(name,S/'competition-original.py');module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
a=read_aam_checkpoint(PACKAGE/'inputs/case64.pkl.gz')
old=module.compete_fragments(a);new=compete_fragments(a,CompetitionConfig())
assert old.counts==new.counts and len(old.repairs)==len(new.repairs) and old.pending==new.pending
for x,y in zip(old.offers,new.offers):
 y=dict(y);y.pop('proposal_mode');assert x==y
for x,y in zip(old.repairs,new.repairs):assert x.graph.to_record()==y.graph.to_record()
print('Local control exactly reproduces production: counts, offers and all repair graphs agree.')
