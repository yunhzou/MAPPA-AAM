"""Produce reviewable complete-study tables; never publish partial numbers."""
from pathlib import Path
import csv,json,hashlib,importlib.metadata,sys,gzip
import run
from aggregate import aggregate, stats
S=Path(__file__).resolve().parent
r=aggregate()
# Save diagnostics even when task completion or hardware consistency is missing.
if r['status']!='all_tasks_processed' or not r['same_cpu_model'] or not r['all_threads_one']:
 run.save(S/'finalization-blocked.json',dict(status='blocked',processed_cases=r['processed_cases'],cpu_models=r['cpu_models'],threads_one=r['all_threads_one'],reason='Complete same-hardware records required before producing final comparison tables.'))
 raise SystemExit(2)
rows=[]
for key,m in r['methods'].items():
 cfg=m['config'];paired=m['paired_completed_api_cpu'];all_api=m['all_attempted_api_cpu'];worker=m['all_attempted_worker_cpu'];o=m['outcomes']
 rows.append(dict(configuration=key,branch_cap=cfg.get('cap','NA'),seed_orders=cfg.get('seed','NA'),sweep=cfg['sweep'],reactions=1851,recovered=o.get('recovered',0),not_recovered=o.get('not_recovered',0),unknown=o.get('unknown',0),recovery_percent=100*o.get('recovered',0)/1851,search_completed_cases=m['completed_search_cases'],common_completed_cases=len(r['common_completed_cases']),paired_mean_cpu=paired['mean_cpu_seconds'] if paired else None,paired_median_cpu=paired['median_cpu_seconds'] if paired else None,paired_p95_cpu=paired['p95_cpu_seconds'] if paired else None,all_attempted_api_n=all_api['n'] if all_api else 0,all_attempted_api_mean_cpu=all_api['mean_cpu_seconds'] if all_api else None,all_attempted_worker_n=worker['n'] if worker else 0,all_attempted_worker_mean_cpu=worker['mean_cpu_seconds'] if worker else None))
with (S/'comparison.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
# Size-oriented policies use explicit atoms, with stored-order tie breaking.
sizes={}
for case in range(1851):
 raw=run.read(run.BASE/'golden-inputs'/str(case)/'input.json')
 sizes[case]=(len(raw['reactant']['elements']),len(raw['product']['elements']))
directional=[]
common_set=set(r['common_completed_cases'])
for key,m in r['methods'].items():
 for policy in ['smaller_first','larger_first','bidirectional']:
  outcomes={'recovered':0,'not_recovered':0,'unknown':0};cpu=[]
  for row in m['per_case']:
   case=row['case'];nr,np=sizes[case];small='R_to_P' if nr<=np else 'P_to_R';large='P_to_R' if small=='R_to_P' else 'R_to_P'
   directions=[small,large] if policy=='bidirectional' else [small if policy=='smaller_first' else large]
   verdicts=[row['directions'][d]['reference_recovery'] for d in directions]
   verdict='recovered' if 'recovered' in verdicts else 'not_recovered' if all(v=='not_recovered' for v in verdicts) else 'unknown'
   outcomes[verdict]+=1
   if case in common_set:cpu.append(sum(row['directions'][d]['api_cpu_seconds'] for d in directions))
  summary=stats(cpu)
  directional.append(dict(configuration=key,policy=policy,reactions=1851,**outcomes,recovery_percent=100*outcomes['recovered']/1851,common_completed_cases=len(cpu),mean_cpu=summary['mean_cpu_seconds'] if summary else None,median_cpu=summary['median_cpu_seconds'] if summary else None,p95_cpu=summary['p95_cpu_seconds'] if summary else None))
with (S/'directional.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(directional[0]));writer.writeheader();writer.writerows(directional)

lines=['# Controlled Golden comparison','',f"CPU model: {next(iter(r['cpu_models']))}. One pinned thread per search. All recovery percentages use 1,851 reactions. The paired completed-search timing cohort contains {len(r['common_completed_cases'])} identical reactions for all configurations.",'','Search workflow CPU includes graph preparation, mapping and native output persistence, and excludes reference verification and initial process imports. MAPPA retains raw compressed graphs; SLAP retains native label information and explicit candidates. Separate all-attempted CPU means include resource-limited calls and describe spent computation, not completed-search latency. See CSV for denominators and completed-case counts.','', '| Configuration | Recovered | Unknown | Paired mean CPU s/reaction | All-attempted mean CPU s/reaction |','|---|---:|---:|---:|---:|']
for row in rows:
 def fmt(x):return 'unavailable' if x is None else f'{x:.3f}'
 lines.append(f"| {row['configuration']} | {row['recovered']}/1851 ({row['recovery_percent']:.2f}%) | {row['unknown']} | {fmt(row['paired_mean_cpu'])} | {fmt(row['all_attempted_api_mean_cpu'])} (n={row['all_attempted_api_n']}) |")
(S/'COMPARISON.md').write_text('\n'.join(lines)+'\n')
packages={}
for name in ['numpy','scipy','sympy','pynauty','z3-solver','rdkit','psutil']:
 try:packages[name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError:packages[name]=None
run.save(S/'runtime-proof.json',dict(python=sys.version,packages=packages,native_binaries={str(p.relative_to(S)):run.sha(p) for p in (S/'engine/src/mappa').glob('*.so')},cpu_models=r['cpu_models']))
run.save(S/'finalization.json',dict(status='tables_ready_for_review',processed_cases=1851,common_completed_cases=len(r['common_completed_cases']),source_manifest_sha256=run.sha(S/'manifest.json'),files={n:run.sha(S/n) for n in ['controlled-results.json.gz','comparison.csv','directional.csv','COMPARISON.md','runtime-proof.json']},scope='These are measured results awaiting manuscript integration and visual review. No paper or Git publication is automatic.'))
print('Comparison tables ready for review',flush=True)
