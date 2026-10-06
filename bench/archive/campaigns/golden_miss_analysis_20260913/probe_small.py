"""Reference-blind sensitivity probes, excluded from benchmark scores."""
import sys,json,time,resource,hashlib
from dataclasses import asdict
from pathlib import Path
ROOT=Path(sys.argv[1]);OUT=Path(sys.argv[2]);CASE=int(sys.argv[3]);DIRECTION=sys.argv[4];VARIANT=sys.argv[5]
sys.path[:0]=[str(ROOT/'engine/src'),str(ROOT/'engine/bench')]
from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,search_aam
from mappa.search_orientation import AAMSearchPlan
from golden_evaluation import evaluate_planned
raw=json.loads((ROOT/f'golden-inputs/{CASE}/input.json').read_text())
problem=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw.get('name',''))
config=AAMSearchConfig(seed_count=10,branch_limit=1000 if VARIANT=='cap1000' else 100,iso_tolerance=.5 if VARIANT=='tol05' else 1.)
reverse=DIRECTION=='P_to_R';plan=AAMSearchPlan(problem,AAMProblem(problem.product,problem.reactant,problem.name) if reverse else problem,config,reverse)
folder=OUT/f'case{CASE}-{DIRECTION}-{VARIANT}';folder.mkdir(parents=True,exist_ok=True)
start=time.perf_counter();cpu=time.process_time()
aam=search_aam(plan.problem,config,execution='reused_native',workers=1,intermediate_dir=folder/'cuts',archive_format='checkpoint')
search_cpu=time.process_time()-cpu
# Reference is loaded only after search; it cannot guide the mapping.
reference=json.loads((ROOT/f'golden-inputs/{CASE}/reference.json').read_text())
evaluation=evaluate_planned(aam,plan,reference['features'],reference['mapping'],seconds=120,query_timeout_ms=5000)
result=dict(case=CASE,direction=DIRECTION,variant=VARIANT,config=asdict(config),metrics=asdict(aam.metrics),evaluation=evaluation,search_cpu_seconds=search_cpu,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,reference_directed=False,included_in_benchmark=False,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(folder/'result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(case=CASE,direction=DIRECTION,variant=VARIANT,verdict=evaluation['reference_recovery'],capped=aam.graph.capped,seconds=result['wall_seconds'])),flush=True)
