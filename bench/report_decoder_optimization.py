"""Summarize/plot an exactness-checked saved-search decoder timing run."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mappa-decoder-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--run',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
repo=Path(__file__).resolve().parents[1]
read=lambda p:json.loads(p.read_text())
old=read(repo/'reports/final_end_to_end_20260916/summary.json')
rows=[read(p) for p in args.run.glob('case*/result.json')]
fresh = bool(rows and 'stages' in rows[0])
if fresh:
    old_support = {}
    with gzip.open(repo/'reports/final_end_to_end_20260916/family-journals.jsonl.gz','rt') as f:
        for line in f:
            record=json.loads(line)
            if record['complete']:
                old_support.setdefault(record['case'], {})[record['family']]=sorted(record['pattern_ids'])
    converted=[]
    for r in rows:
        journal=[json.loads(line) for line in (args.run/f"case{r['case']}/families.jsonl").read_text().splitlines()]
        checked=old_support.get(r['case'],{})
        converted.append(dict(case=r['case'],complete=r['complete'],decode=r['stages']['decode'],
            families=r['families'],pattern_ids=sorted(r['observed_classes']),
            same_class_ids=r['validation']['same_class_ids'],same_family_count=r['validation']['same_family_count'],
            compared_family_supports=len(checked),same_compared_family_supports=all(
                journal[i]['complete'] and sorted(journal[i]['pattern_ids'])==ids for i,ids in checked.items())))
    rows=converted
rows.sort(key=lambda r:r['case'])
assert len(rows)==140 and all(all(r[k] for k in ('complete','same_class_ids','same_family_count',
    'same_compared_family_supports')) for r in rows)
assert sum(len(r['pattern_ids']) for r in rows)==300
previous={r['execution']['case']:r['result']['stages']['decode']['cpu_seconds']
          for r in old['rows'] if r['result'] and r['execution']['complete']}
new={r['case']:r['decode']['cpu_seconds'] for r in rows}
paired=sorted(previous)
def stats(v):
    return dict(mean=float(np.mean(v)),median=float(np.median(v)),p95=float(np.percentile(v,95)),
                total=float(np.sum(v)),maximum=float(np.max(v)))
summary=dict(fresh_end_to_end=fresh,completed=140,all_class_sets_match=True,
    compared_family_supports=sum(r['compared_family_supports'] for r in rows),
    all_compared_family_supports_match=True,all_families=sum(r['families'] for r in rows),
    paired_cases=paired,before_paired=stats([previous[c] for c in paired]),
    after_paired=stats([new[c] for c in paired]),after_all=stats(list(new.values())),
    manifest=read(args.run/'manifest.json'),campaign=read(args.run/'campaign.json'),rows=rows)
args.output.mkdir(parents=True,exist_ok=True)
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with gzip.open(args.output/'family-journals.jsonl.gz','wt') as out:
    for r in rows:
        for line in (args.run/f"case{r['case']}/families.jsonl").read_text().splitlines():
            out.write(json.dumps(dict(case=r['case'],**json.loads(line)),separators=(',',':'))+'\n')
a,b=summary['before_paired'],summary['after_paired']
text=f'''# MAPPA decoder optimization

All 140 saved searches decode completely with the same 300 reaction-local event classes.
All {summary['compared_family_supports']:,} previously completed per-family class-support records match.
The two previously interrupted runs now finish. The saved catalogue still has {summary['all_families']:,} families.

| CPU seconds | Before (same 138 cases) | After (same 138 cases) | After (all 140) |
|---|---:|---:|---:|
'''
for key,label in [('mean','Mean'),('median','Median'),('p95','95th percentile'),('total','Total')]:
    text+=f"| {label} | {a[key]:.3f} | {b[key]:.3f} | {summary['after_all'][key]:.3f} |\n"
text+=f'''
New timing source: {'a fresh search-plus-decoding run; full stage timings are preserved in the accompanying end-to-end report' if fresh else 'saved-search decoder replay'}.

Paired decoder CPU reduction: {100*(1-b['total']/a['total']):.1f}% ({a['total']/b['total']:.2f}x faster in aggregate).
Case 25 now takes {new[25]:.3f} decoder CPU seconds; case 31 takes {new[31]:.3f}.
Their old attempts hit the 300-second **whole-process** watchdog; these are censored outcomes, not completed decoder times.

Four independent reaction processes, one numerical thread each; decoder CPU is not divided by four.
Both decoder measurements include lazy initialization and per-family journaling; search and catalogue construction are excluded from decoder stage times. The summary records whether the new decoding followed a fresh search or a saved checkpoint.
Thresholds (.5/.3), matching tolerance, saved families, event windows, and exact class identity are unchanged. No candidate cap or bijection/group-element enumeration was added.

Changes:
- Reuse exact raw-weight pair response tables and symbolic subset selectors.
- Express event totals as Boolean cardinality constraints instead of arithmetic search.
- Omit a redundant global injectivity constraint: each pool/group action is already a permutation.
- Share source edges across families; screen absent bonds in NumPy before scalar witness validation.

The 16,384-entry local-table cache is a memory policy, not a solution limit: evicted entries are recomputed exactly.
All original families and symmetry-query information remain available.

Reproduction:

```sh
python bench/benchmark_saved_decoding.py --archives /path/to/final_end_to_end_run --output /path/to/new_decode_run --workers 4
python bench/report_decoder_optimization.py --run /path/to/new_decode_run --output reports/decoder_optimization_20260916
```

Source hashes, per-reaction class IDs, timings, and validation flags are in `summary.json`; exact per-family results are in `family-journals.jsonl.gz`.
'''
(args.output/'README.md').write_text(text)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,(ax,bx)=plt.subplots(1,2,figsize=(11.8,5.0),gridspec_kw={'width_ratios':[1.45,1]})
fig.subplots_adjust(left=.085,right=.975,top=.77,bottom=.23,wspace=.37)
fig.text(.085,.93,'MAPPA decoding: less repeated work, same candidates',fontsize=19,weight='bold',color='#183747')
fig.text(.085,.855,f"140/140 complete · 300 event classes preserved · {summary['compared_family_supports']:,} family-support checks agree",fontsize=11,color='#506772')
x=np.array([previous[c] for c in paired]);y=np.array([new[c] for c in paired])
ax.scatter(x,y,s=25,color='#20858B',alpha=.75)
lo=min(x.min(),y.min())*.7;hi=max(x.max(),y.max())*1.4
ax.plot([lo,hi],[lo,hi],color='#9AAAB0',ls='--',lw=1)
ax.set(xscale='log',yscale='log',xlim=(lo,hi),ylim=(lo,hi),xlabel='Previous decoder CPU seconds',ylabel='Optimized decoder CPU seconds',title='Same 138 completed reactions')
ax.text(.04,.91,'Below the line = faster',transform=ax.transAxes,color='#20858B',fontsize=10)
ax.grid(alpha=.12)
bx.bar([0,1],[a['mean'],b['mean']],color=['#A6AFBA','#20858B'],width=.55)
for i,val in enumerate([a['mean'],b['mean']]):bx.text(i,val+.12,f'{val:.2f} s',ha='center',weight='bold',fontsize=14)
bx.set(xticks=[0,1],xticklabels=['Previous','Optimized'],ylabel='Mean decoder CPU seconds',ylim=(0,a['mean']*1.3),title=f"{a['total']/b['total']:.2f}× aggregate speedup")
bx.grid(axis='y',alpha=.12);bx.set_axisbelow(True)
fig.text(.085,.10,f"Previously timed out: case 25 now {new[25]:.1f} s; case 31 now {new[31]:.1f} s (decoder CPU).",fontsize=11)
fig.text(.085,.047,('Fresh search + decoding' if fresh else 'Saved searches')+' · cap 2,000 · one-seed sweep · competition off · identical event windows · four serial workers',fontsize=9,color='#506772')
for ext in ('png','svg'):fig.savefig(args.output/f'decoder-speedup.{ext}',dpi=180,facecolor='white')
print(json.dumps({k:v for k,v in summary.items() if k not in ('manifest','campaign','rows','paired_cases')},indent=2))
