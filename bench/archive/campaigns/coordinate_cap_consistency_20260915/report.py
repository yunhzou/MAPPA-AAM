"""Write a transparent cap comparison from the saved paired evidence."""
from pathlib import Path
import json,collections
from analyze import analyze
S=Path(__file__).resolve().parent
r=analyze();data=json.loads((S/'comparison.json').read_text());controls=json.loads((S/'controls.json').read_text())['per_case']
def cases(v):return ', '.join(map(str,v)) if v else 'none'
def changes(v):return '; '.join(f"case {x['case']}: {x['count']}" for x in v) if v else 'none'
lines=['# Branch-cap comparison on 140 coordinate reactions','',('Complete paired evaluation.' if r['processed']==140 and not r['unresolved'] else f"Interim result: {r['processed']}/140 processed; {r['complete_pairs']} paired decoding windows complete."),'','Both configurations use one seed, the cut sweep, fragment competition, and the same event windows. Only the branch cap changes: 100 or 2,000. The cap applies to initial growth and competition. Search and event decoding remain separate.','', '| Result | Cap 100 | Cap 2,000 |','|---|---:|---:|']
for label,key in [('Complete saved-family decoding windows','complete_saved_full_window'),('Reactions with a certified minimum in the saved full families','minimum_proven'),('Reactions without a complete mapping','no_full_families')]:
 vals=[sum(bool(row['caps'][str(cap)]['decoded'] and row['caps'][str(cap)]['decoded'][key]) for row in data['per_case']) for cap in (100,2000)]
 lines.append(f'| {label} | {vals[0]} | {vals[1]} |')
lines+=['',f"Equal minimum counts on {r['same_minimum']} reactions. Higher cap-100 minima: {cases(r['higher_minimum'])}; lower: {cases(r['lower_minimum'])}. Cap-100 searches without complete mappings: {cases(r['no_full_mapping']['100'])}. Incomplete paired decoding: {cases(r['unresolved'])}.",'',f"Cap-2,000 minimum patterns lost at cap 100: {changes(r['minimum_patterns_lost'])}. All-window pattern losses: {changes(r['window_patterns_lost'])}. Additional cap-100 window patterns: {changes(r['window_patterns_added'])}.",'',f"Fresh cap-2,000 cases disagreeing with saved minimum counts or omitting a saved minimum pattern: {cases(r['fresh2000_disagrees_with_saved'])}. Full-window differences from the saved cap-2,000 output: {cases([x['case'] for x in r['fresh2000_window_differences']])}.",'','## Comparison with the saved SLAP results','','Mappings in this collection are unverified. These results compare event counts and patterns, not mapping accuracy. Minima refer to returned candidates, not global optima.','', '| MAPPA cap | Comparator | Fewer events | Equal | More | No full mapping | Unresolved | Comparator minimum patterns covered in window |','|---|---|---:|---:|---:|---:|---:|---:|']
for cap in ('100','2000'):
 for method in ('native_slap','slap_sweep'):
  counts=collections.Counter();covered=total=0
  for row in data['per_case']:
   d=row['caps'][cap]['decoded'];c=controls[row['case']]['comparators'][method]
   total+=len(c['ids']);covered+=len(set(c['ids'])&set(d['patterns'] if d else []))
   if not d or not d['complete_saved_full_window']:counts['unresolved']+=1
   elif d['no_full_families']:counts['no_mapping']+=1
   elif not d['minimum_proven']:counts['unresolved']+=1
   else:counts['lower' if d['minimum']<c['minimum'] else 'higher' if d['minimum']>c['minimum'] else 'equal']+=1
  lines.append(f"| {cap} | {method} | {counts['lower']} | {counts['equal']} | {counts['higher']} | {counts['no_mapping']} | {counts['unresolved']} | {covered}/{total} |")
lines+=['','## Timing','',f"CPU: {', '.join(data['cpu_models'])}. {r['timing']['scope']}",'','| Stage | Cap 100 mean CPU s/reaction | Cap 2,000 mean CPU s/reaction |','|---|---:|---:|']
for phase in ('search_internal','competition_internal','decode','total_worker'):
 a,b=[r['timing'][str(c)][phase] for c in (100,2000)]
 if a['n'] and b['n']:lines.append(f"| {dict(search_internal='Sweep search',competition_internal='Fragment competition',decode='Decoding worker',total_worker='Total worker CPU')[phase]} (n={a['n']}) | {a['mean']:.3f} | {b['mean']:.3f} |")
lines+=['', 'All processed worker CPU (including empty results and all recorded decoder passes): '+ '; '.join(f"cap {cap}: {r['all_processed_worker_cpu'][cap]['total']:.3f} seconds over {r['all_processed_worker_cpu'][cap]['n']} reactions" for cap in ['100','2000'])+'.']
lines+=['','The worker total includes process imports, so it is not the sum of import-excluded search/competition and decoder rows. Timing uses identical reactions with completed mappings at both caps; failed searches and unfinished decoding are excluded from this paired latency summary. Raw stage records preserve their computational cost. SLAP times from a different CPU are not used to construct a speed ratio.','','## Interpretation','','Case 123 reaches the cap before obtaining a complete parent mapping at cap 100. Competition has no complete parent to revise, so it cannot repair this search. At cap 2,000 the initial search returns 73 complete terminals and the final pipeline recovers the saved three-event pattern. This is a search-budget failure, not an unresolved decoder or proof that the reaction has no mapping.','','Decoding certificates cover the retained complete families within the common event window. Neither branch cap establishes exhaustive search. A lower cap may also change which competition candidates receive its fixed budget, so the decoded sets are compared directly rather than assumed to be nested.','','Case 125 is a second confirmed empty search at cap 100; cap 2,000 returns its saved four-event solution. Both cases must count as search failures under a uniform cap-100 protocol.','','Source and input hashes are frozen in manifest.json. The experimental driver is run.py; analyze.py interprets complete empty searches separately from unresolved decoding. Raw mappings and private execution details are excluded from this report.']
(S/'RESULTS.md').write_text('\n'.join(lines)+'\n')
