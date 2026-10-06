"""Small guarded integration check against existing one-seed cut archives."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import resource
import sys
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
sys.path[:0]=[str(REPO/'src'),str(REPO/'bench')]
from golden_checkpoint_evaluation import evaluate_checkpoints
from mappa import AAMProblem, MolecularEndpoint, AAMSearchConfig, AAMSearchPlan

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('campaign',type=Path);parser.add_argument('case',type=int);parser.add_argument('direction',choices=['R_to_P','P_to_R']);parser.add_argument('output',type=Path)
 args=parser.parse_args()
 raw=json.loads((args.campaign/f'golden-inputs/{args.case}/input.json').read_text())
 ref=json.loads((args.campaign/f'golden-inputs/{args.case}/reference.json').read_text())
 problem=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw.get('name',''))
 cfg=AAMSearchConfig(seed_count=1,branch_limit=100)
 reverse=args.direction=='P_to_R'
 plan=AAMSearchPlan(problem,AAMProblem(problem.product,problem.reactant,problem.name) if reverse else problem,cfg,reverse)
 result=evaluate_checkpoints(args.campaign/f'runs/golden/seed1/case{args.case}/{args.direction}/cuts',plan,ref['features'],ref['mapping'])
 result.update(case=args.case,direction=args.direction,seed_count=1,config=asdict(cfg),peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024))
 assert result['reference_recovery']=='not_recovered', result
 args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ('case','direction','reference_recovery','wall_seconds','cpu_seconds','peak_mib')}),flush=True)
