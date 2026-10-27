"""Capture one peak-fragment witness, checking the saved capped growth outcome."""
import json,tarfile,tempfile,sys,hashlib
from pathlib import Path
from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,aam
from mappa.artifacts import read_raw_cut
from mappa.growth import native

archive,output,seed=sys.argv[1:];seed=int(seed)
with tarfile.open(archive) as t, tempfile.TemporaryDirectory(dir='/tmp') as td:
 manifest=json.load(t.extractfile('archive-manifest.json'))['files']
 def read(name):
  data=t.extractfile(name).read();assert hashlib.sha256(data).hexdigest()==manifest[name]['sha256'];return data
 raw=json.loads(read('input.json')); task=json.loads(read('task.json'))
 path=Path(td)/'cut.raw.pkl.gz';path.write_bytes(read('cuts/cut_00000.raw.pkl.gz'));saved=read_raw_cut(path)
context=saved.contexts[seed]
endpoints=[MolecularEndpoint(**raw[k]) for k in ('reactant','product')]
if task['direction']=='P_to_R':endpoints.reverse()
config=AAMSearchConfig(**task['config']);problem=AAMProblem(*endpoints,raw['name'])
aam._initialize_search(problem,config,'reused_native');view=aam._SEARCH_REPAIR.for_cut(context.cuts)
r=native.source_graph(view.source);p=native.target_graph(view.target,aam._SEARCH_CONTEXT[3])
args=(r.graph,p.graph,context.seed_order[0],[-1]*len(r.nodes),config.graph_floor,config.iso_tolerance,1,config.branch_limit,[],[],False,config.random_seed)
plain=native._engine.grow_island(*args)
observed=native._engine.grow_island(*args,observe=True)
snapshot=observed.pop('observation');assert plain==observed
stop=next(s for s in saved.stops if s.seed==context.seed_order[0])
assert plain['capped'] and plain['cap_count']==stop.count
assert len(saved.contexts)==10
result=dict(snapshot=snapshot,normal=plain,seed_index=seed,seed_atom=context.seed_order[0],seed_order=context.seed_order,direction=task['direction'],config=task['config'],observation_equality=True,saved_stop_equal=True)
Path(output).write_text(json.dumps(result)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('snapshot','seed_order','config')}),flush=True)
print('peak',snapshot['candidate_count'],'fragment atoms',len(snapshot['fragment']),flush=True)
