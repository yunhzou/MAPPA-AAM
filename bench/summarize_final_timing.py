"""Summarize measured stages without treating watchdogs as completed runtimes."""
from pathlib import Path
import hashlib,json,gzip,statistics
import numpy as np
import argparse
parser=argparse.ArgumentParser(description='Summarize fresh final MAPPA end-to-end timings, retaining censored attempts.')
parser.add_argument('run',type=Path)
args=parser.parse_args()
W=args.run
def read(p): return json.loads(p.read_text())
def stats(xs):
 return dict(n=len(xs),mean=float(np.mean(xs)),median=float(np.median(xs)),p95=float(np.percentile(xs,95)),maximum=max(xs),total=sum(xs)) if xs else None
rows=[]
for p in sorted(W.glob('case*/execution.json'),key=lambda p:int(p.parent.name[4:])):
 e=read(p);rpath=p.parent/'result.json';ppath=p.parent/'progress.json'
 row=dict(execution=e,result=read(rpath) if rpath.exists() else None,progress=None if rpath.exists() or not ppath.exists() else read(ppath))
 journal=p.parent/'families.jsonl'
 if journal.exists():
  records=[]
  for line in journal.read_text().splitlines():
   try: records.append(json.loads(line))
   except json.JSONDecodeError: break
  row['journal_summary']=dict(processed=len(records),certified=sum(x['complete'] for x in records),classes=sorted({i for x in records for i in x['pattern_ids']}))
 rows.append(row)
complete=[r['result'] for r in rows if r['execution']['complete']]
for row in rows:
 if row['execution']['complete']:
  r=row['result']
  assert r['complete'] and r['processed_families']==r['completed_families']==r['families']
  assert row['journal_summary']['certified']==r['families']
  assert set(row['journal_summary']['classes'])==set(r['observed_classes'])
  for unit in ['cpu_seconds','wall_seconds']:
   assert sum(d[unit] for d in r['stages'].values()) <= r['end_to_end'][unit]+1e-6
summary=dict(attempted=len(rows),requested=140,completed=len(complete),complete_cases=[r['case'] for r in complete],
 incomplete_cases=[r['execution']['case'] for r in rows if not r['execution']['complete']],
 errors=[r['execution'] for r in rows if r['execution']['status']=='error'],
 all_completed_match=all(all(r['validation'][k] for k in ['same_class_ids','same_family_count']) for r in complete),
 complete_cohort={k:{unit:stats([r['end_to_end' if k=='end_to_end' else 'stages'][unit] if k=='end_to_end' else r['stages'][k][unit] for r in complete]) for unit in ['cpu_seconds','wall_seconds']} for k in ['search','catalogue','decode','serialize','end_to_end']},
 measured_search_cohort={unit:stats([d['stages']['search'][unit] for row in rows if (d:=row['result'] or row['progress']) and 'search' in d['stages']]) for unit in ['cpu_seconds','wall_seconds']},
 process_wall_sum_seconds=sum(r['execution']['process_wall_seconds'] for r in rows),
 sampled_process_cpu_sum_seconds=sum(r['execution']['sampled_process_cpu_seconds'] for r in rows),
 peak_worker_mib=max((r['execution']['peak_mib'] for r in rows),default=0),
 scope='Completed-cohort distributions include fresh input loading, search/checkpoint I/O, deduplication, public event decoding, result serialization, and light journal overhead; exclude setup/imports and class-ID comparison. Incomplete attempts are censored and excluded from completed-latency averages, listed separately.',
 rows=rows)
summary['manifest']=read(W/'manifest.json')
summary['campaign']=read(W/'campaign.json')
(W/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# Final MAPPA end-to-end timing', '', f"Fresh attempts: {len(rows)}/140; completed within the five-minute watchdog: {len(complete)}.", '', 'One seed, single-edge cut sweep, cap 2,000, R-to-P, competition off. Event windows match the published 140-case comparison. No candidate-count cap. Four processes, each using one numerical thread.', '', '**The following distributions cover only completed reactions. They are not full-dataset completed-runtime estimates when any reaction is interrupted.**', '', '| Stage | Mean CPU s | Median CPU s | 95th percentile CPU s | Median wall s |', '|---|---:|---:|---:|---:|']
for key,label in [('search','Search'),('catalogue','Final-family catalogue'),('decode','Event decoding'),('serialize','Output serialization'),('end_to_end','End to end')]:
 d=summary['complete_cohort'][key];c=d['cpu_seconds'];w=d['wall_seconds']
 if c:lines.append(f"| {label} | {c['mean']:.3f} | {c['median']:.3f} | {c['p95']:.3f} | {w['median']:.3f} |")
lines += ['', 'End to end includes input loading, search/checkpoint I/O, catalogue construction, decoding, candidate serialization, and light progress journaling. Initial imports and reference-set validation are excluded. The decoder includes lazy dependency initialization. Stage medians and percentiles are not additive.', '', f"Incomplete cases: {summary['incomplete_cases']}. See summary.json for resource-stop status, incurred process time, partial progress, and completed family certificates. An interrupted attempt is not a completed 300-second result.", '', f"All completed class sets and family counts match the archived baseline: {summary['all_completed_match']}. No saved search or decoder outputs were reused in the timed pipeline.", '', 'Reproduction (from the repository root):', '', '```sh', 'python native/build_engine.py', 'python bench/final_end_to_end.py --inputs /path/to/full140_inputs --output /path/to/new_run --cases all --workers 4', 'python bench/summarize_final_timing.py /path/to/new_run', '```', '', 'The run manifest records source, native-binary and input hashes. Raw per-case search checkpoints, candidates and family journals remain in the run directory. This report preserves per-case timing and verification summaries.']
# Preserve full-attempt spending alongside the completed-cohort distributions.
spent=summary['sampled_process_cpu_sum_seconds']
insert=lines.index('**The following distributions cover only completed reactions. They are not full-dataset completed-runtime estimates when any reaction is interrupted.**')
lines[insert:insert]=['## All attempted work, including timeouts', '',
 f'Observed process CPU spent: approximately {spent/60:.2f} CPU minutes, or {spent/len(rows):.2f} CPU seconds per attempted reaction. This sampled process scope includes startup and interrupted calls. It is not the cost of completing all reactions.', '']
if not summary['incomplete_cases']:
 lines[insert] = '## Full campaign cost'
 lines[insert+2] = (f'Observed process CPU spent: approximately {spent/60:.2f} CPU minutes, '
  f'or {spent/len(rows):.2f} CPU seconds per reaction. All reactions completed. '
  'This sampled process scope includes startup; the pipeline timings below exclude initial imports.')
 lines = [('**The following distributions cover all 140 reactions.**' if line.startswith('**The following distributions cover only') else
           'All 140 reactions completed; there were no watchdog or resource stops.' if line.startswith('Incomplete cases: []') else line)
          for line in lines]
(W/'README.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['rows','complete_cohort','complete_cases','manifest','campaign']},indent=2))
for k,v in summary['complete_cohort'].items(): print(k,v)
