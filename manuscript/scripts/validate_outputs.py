"""Validate the final manuscript's numerical sources and compiled references."""
from pathlib import Path
import hashlib,json,re
from publication_privacy import check_publication
from pypdf import PdfReader
MAN=Path(__file__).resolve().parents[1];E=MAN/'evidence'
check_publication('manuscript.pdf',(MAN/'manuscript.pdf').read_bytes())
def read(n):return json.loads((E/n).read_text())
from build_molecule_figure import check_example
assert read('molecule_example.json') == json.loads(json.dumps(check_example()))
from build_alternatives_figure import check_example as check_alternatives
check_alternatives(MAN)
from build_alternatives_figure import check_case
check_case(MAN)
uncut=read('unswept.json')
assert uncut['denominator']==1851 and uncut['seed_count']==1
for name in ['mappa','slap']:
 d=uncut['methods'][name]
 assert sum(d['counts'].values())==1851
 assert d['counts']['recovered']==len(set(d['recovered_cases']))
 assert d['sweep_recovered']==len(d['recovered_cases'])+len(d['added_by_sweep'])
competitors=read('competitors.json')
assert competitors['denominator']==1851 and competitors['rescored_with_current_evaluator'] and competitors['unchanged_case_outcomes']
assert len(competitors['methods'])==7
for d in competitors['methods']:
 assert d['total']==1851 and sum(d['statuses'].values())==1851
 assert d['first_correct']==len(set(d['first_correct_cases']))
 assert d['any_correct']==len(set(d['any_correct_cases']))
 assert set(d['first_correct_cases'])<=set(d['any_correct_cases'])<=set(range(1851))
 assert len(d['failed_cases'])==sum(n for k,n in d['statuses'].items() if k!='mapped')
 assert d['invalid_candidates']==sum(d['invalid_reasons'].values())
seed=read('seed_comparison.json');slap=read('slap_sweep.json');flat=read('final_dedup.json')
misses=read('golden_miss_analysis.json')
assert misses['nonrecovered']==11 and len(misses['rows'])==11
assert {r['case'] for r in misses['rows']}=={r['case'] for r in seed['methods']['seeds10']['per_case'] if r['outcome']=='not_recovered'}
assert sum(r['classification']=='extra_matched_pair' for r in misses['rows'])==4
assert all(r['retained_pair_counts']==[r['reference_pairs']+1] for r in misses['rows'] if r['classification']=='extra_matched_pair')
assert misses['all_reference_conversions_verified'] and misses['cases_recovered_by_any_comparator']==[871]
assert all(t['native_python_graph_equal'] for t in misses['causal_traces'])

assert seed['fresh'] and slap['fresh']
keys=['seeds1','seeds2','seeds3','seeds10']
assert all(seed['methods'][k]['golden_cases']==1851 for k in keys)
assert all(sum(seed['methods'][k]['golden_outcomes'].values())==1851 for k in keys)
for k in keys:
 d=seed['methods'][k]
 assert len(d['per_case'])==1851 and {r['case'] for r in d['per_case']}==set(range(1851))
 assert all(r['outcome'] in {'recovered','not_recovered','unknown'} for r in d['per_case'])
 for status in ('recovered','not_recovered','unknown'):
  assert d['golden_outcomes'][status]==sum(r['outcome']==status for r in d['per_case'])
 if d['common_mean_cpu_seconds'] is not None:
  assert abs(d['common_mean_cpu_seconds']-sum(d['per_case'][c]['search_cpu_including_io'] for c in seed['common_case_indices'])/len(seed['common_case_indices']))<1e-8
  assert all(d['per_case'][c]['search_complete'] and d['per_case'][c]['search_hosts']==['Mac'] for c in seed['common_case_indices'])
assert sum(slap['outcomes'].values())==1851
assert slap['sweep_union_recovered']==slap['outcomes']['recovered']
assert flat['competition'] is False
assert len(flat['per_case'])==140 and all(r['complete'] for r in flat['per_case'])
assert flat['flat_families']==120052==sum(r['flat_saved_families'] for r in flat['per_case'])
assert flat['event_classes']==300==sum(len(r['class_ids']) for r in flat['per_case'])
assert sum(flat['window_distribution'].values())==140
sources=json.loads((E/'paper_sources.json').read_text())
assert sources['fresh_campaign_complete']
for row in sources['snapshots']:
 assert hashlib.sha256((E/row['snapshot']).read_bytes()).hexdigest()==row['sha256'],row['snapshot']
tex='\n'.join(p.read_text() for p in (MAN/'includes').glob('*.tex') if p.name in ['paper.tex','supplement.tex','include-abstract.tex'])
for obsolete in ['Adaptive no-sweep','Separate experimental versions','Holdout follow-ups','earlier publication engine','N_{\\mathrm{order\\ changed}}']:
 assert obsolete not in tex,obsolete
paper=(MAN/'includes/paper.tex').read_text()
assert paper.index(r'\label{fig:alternatives}') < paper.index(r'\section{Evaluation}')
assert 'Final representation and complete window decoding' not in paper
bib=(MAN/'references.bib').read_text()
for group in re.findall(r'\\cite\w*\{([^}]+)\}',tex):
 for key in group.split(','):assert re.search(r'@\w+\{'+re.escape(key)+',',bib),key
figs=re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',tex)
assert set(figs)=={'figs/fig1_algorithm.pdf','figs/fig2_golden.pdf','figs/fig3_coordinate.pdf','figs/fig4_alternatives.pdf'}
for fig in figs:
 p=MAN/fig
 assert len(PdfReader(p).pages)==1
 for ext in ['svg','png']:assert p.with_suffix('.'+ext).is_file()
log=(MAN/'build/preprint.log').read_text()
assert not re.search(r'(?:Citation|Reference).*undefined|There were undefined|Overfull \\[hv]box|Missing character:',log), 'Inspect TeX log'
r=PdfReader(MAN/'manuscript.pdf');pages=[p.extract_text() for p in r.pages];alltext='\n'.join(pages)
assert '??' not in alltext
assert not re.search(r'\bMac\b|same-host|worker count|ran concurrently',alltext), 'Keep execution details in the reproducibility records'
for value in [f"{seed['methods']['seeds1']['golden_recovery_percent']:.2f}",f"{seed['methods']['seeds3']['golden_recovery_percent']:.2f}",'120,052','300','162','168']:
 assert value in alltext,value
assert 'Golden reference recovery and runtime' in alltext
assert 'First correct' not in alltext and 'two ten-order cases remain unresolved' not in alltext
comparison=(MAN/'includes/generated-competitor-table.tex').read_text()
assert comparison.count('MAPPA &')==5
assert comparison.count('SLAP,')==3
assert comparison.count('One bijection')==5
assert r'\multicolumn' not in comparison
assert 'SLAP, union' not in comparison and 'SLAP, weighted' not in comparison
assert comparison.count('$^{a}$')==6
assert 'tab:default-times' not in tex
assert not (MAN/'includes/generated-comparator-timing-table.tex').exists()
assert 'generated-seed-table' not in tex
assert not (MAN/'includes/generated-seed-table.tex').exists()
merged=read('golden_cap_ablation.json')
merged_costs={r['key']:r for r in merged['rows']}
merged_calls={r['key']:r for r in read('timing_comparison.json')['default_comparators']}
merged_audit=json.loads((MAN/'build/golden-comparison-table.json').read_text())
assert len(merged_audit['rows'])==13 and merged_audit['no_cross_group_speed_ranking']
for row in merged_audit['rows']:
 expected=(merged_costs[row['timing_key']]['paired_mean'] if row['timing_group']=='A'
           else merged_calls[row['timing_key']]['mean_cpu_seconds'])
 assert row['mean_cpu_seconds']==expected
 assert f"{expected:.3f}" in comparison
 assert row['denominator']==1851
assert 'AMD EPYC 9J14' in alltext and 'CPU model unrecorded' in alltext.replace('\n',' ')
assert '--' not in comparison

for fig in figs:
 assert 'MAPPA' in PdfReader(MAN/fig).pages[0].extract_text(), fig
for value in ['1,489','1,661','80.44','89.74']:assert value in alltext
for d in competitors['methods']:
 if d['method']!='slap_weighted':assert f"{d['any_correct']:,}" in alltext

direction=read('direction_recovery.json')
assert direction['no_search_reruns'] and direction['denominator']==1851
for key,d in direction['methods'].items():
 assert len(d['per_case'])==1851
 for group,g in d['groups'].items():
  rs=[r for r in d['per_case'] if group=='all' or r['cohort']==group]
  assert g['n']==len(rs)
  for policy in ['smaller_first','larger_first','bidirectional']:
   assert g[policy]['cases']==[r['case'] for r in rs if r[policy]=='recovered']
   assert g[policy]['recovered']==len(g[policy]['cases']) and g[policy]['unknown']==0
 assert d['groups']['all']['bidirectional']['recovered']==seed['methods'][key]['golden_outcomes']['recovered']
for value in ['1,803','1,779','97.41','96.11','Smaller-first','Larger-first','(default)']:assert value in alltext,value

timing=read('timing_comparison.json')
assert timing['searches_rerun']==0 and timing['common_case_indices']==seed['common_case_indices']
import statistics
for k,d in timing['methods'].items():
 rows=d['per_case'];assert [r['case'] for r in rows]==seed['common_case_indices']
 xs=sorted(r['cpu_seconds'] for r in rows);st=d['stats'];assert st['n']==1821
 assert abs(st['mean']-statistics.mean(xs))<1e-8 and abs(st['median']-statistics.median(xs))<1e-8
 # 1821 samples give an integer index under the stated linear percentile rule.
 assert abs(st['p95']-xs[1729])<1e-8
for prefix,key in [('mappa1','seeds1'),('mappa2','seeds2')]:
 assert abs(timing['methods'][prefix+'_bidirectional']['stats']['mean']-seed['methods'][key]['common_mean_cpu_seconds'])<1e-8
for v in ['1.141','1.377','2.518','3.329','1.490','4.819','95th pct.']:assert v in alltext,v
for d in timing['default_comparators']:
 source=next(r for r in competitors['methods'] if r['method']==d['key'])
 assert d['calls']==source['statuses']['mapped']
 assert abs(d['mean_wall_seconds']*d['calls']-source['original_mapping_timing']['successful_mapping_wall_sum_seconds'])<1e-8
coord=read('coordinate_minima.json')
assert coord['cases']==140 and coord['mapping_runs']==0
assert coord['mappa_minimum_patterns']==166==sum(len(r['mappa_minimum_ids']) for r in coord['per_case'])
assert coord['comparisons']['slap_sweep']['count_relations']=={'equal':136,'lower':4}
assert coord['comparisons']['native_slap']['count_relations']=={'equal':125,'lower':15}
for method,d in coord['comparisons'].items():
 rs=[r['comparators'][method] for r in coord['per_case']]
 assert d['shared_minimum_patterns']==sum(len(r['shared_minimum_ids']) for r in rs)
 assert d['equal_minimum_patterns']==sum(len(r['ids']) for r in rs if r['count_relation']=='equal')
for d in coord['timing'].values():
 assert len(d['per_case'])==140
 assert abs(d['stats']['mean']-sum(r['cpu_seconds'] for r in d['per_case'])/140)<1e-7
for v in ['158/164','138/140','166','1.047','1.337','0.053']:assert v in alltext,v
assert 'not a mapping-accuracy benchmark' in alltext.replace('\n',' ')
# Prior-method roles and baseline-normalized timing are derived, not new runs.
local=next(d for d in competitors['methods'] if d['method']=='localmapper')
assert local['any_correct']==max(d['any_correct'] for d in competitors['methods'] if not d['method'].startswith('slap'))
assert 'CPU/SLAP' in alltext and 'LocalMapper' in alltext and 'prior' in alltext
for prefix in ['mappa1','mappa2','slap']:
 for policy in ['smaller_first','larger_first','bidirectional']:
  ratio=timing['methods'][prefix+'_'+policy]['stats']['mean']/timing['methods']['slap_'+policy]['stats']['mean']
  assert f'{ratio:.2f}' in (MAN/'includes/generated-direction-table.tex').read_text()
assert min(timing['default_comparators'],key=lambda d:d['mean_wall_seconds'])['key']=='rxnmapper'
assert min(timing['default_comparators'],key=lambda d:d['mean_cpu_seconds'])['key']=='rxnmapper'

# The coordinate cap paragraph carries forward search completion only.
assert '123 and 125' in alltext.replace('\n',' ')
assert not re.search(r'\b(?:competition|takeover|competing)\b',tex,re.I)
assert 'mappa_competition' not in coord['timing']
assert 'mappa_decoding' not in coord['timing']  # no complete baseline timing pass
assert coord['competition'] is False
for row,cert in zip(coord['per_case'],flat['per_case']):
 assert row['case']==cert['case']
 for other in row['comparators'].values():
  assert set(other['covered_in_mappa_window'])==set(other['ids'])&set(cert['class_ids'])

# Same-CPU Golden cap ablation; unresolved is not a verified exclusion.
ablation=read('golden_cap_ablation.json')
assert ablation['cases']==1851 and len(ablation['common_cases'])==1403 and len(ablation['rows'])==12
for row in ablation['rows']:
 assert sum(len(v) for v in row['cases'].values())==1851
 assert len(row['paired_cpu'])==1403 and len(row['selected_attempt_cpu'])==1851
 assert abs(statistics.mean(row['paired_cpu'])-row['paired_mean'])<1e-9
 assert abs(statistics.mean(row['selected_attempt_cpu'])-row['selected_mean'])<1e-9
 assert f"{row['paired_mean']:.3f}" in (MAN/'includes/generated-cap-ablation-table.tex').read_text()
assert all(not r['gains'] and not r['losses'] for r in ablation['paired_cap_audit'])
assert ablation['paired_cap_audit'][-1]['cap100_recovered_now_unknown']==[833,850,1568,1691,1786]
for word in ['1,403','1,835','99.14','43.843','Spent CPU','19 confirmed recoveries']:
 assert word in alltext.replace('\n',' '),word

# Fresh end-to-end timing explicitly retains all interrupted attempts.
e2e=read('end_to_end_timing.json')
assert e2e['attempted']==140 and e2e['completed']+len(e2e['incomplete_cases'])==140
assert e2e['all_completed_match'] and not e2e['errors']
for key,units in e2e['complete_cohort'].items():
 for unit,d in units.items():
  values=[(r['result']['end_to_end'] if key=='end_to_end' else r['result']['stages'][key])[unit] for r in e2e['rows'] if r['execution']['complete']]
  assert d['n']==len(values)==e2e['completed'] and abs(d['mean']-sum(values)/len(values))<1e-8
assert 'MAPPA search plus bond-event decoding' in alltext
assert f"{e2e['complete_cohort']['end_to_end']['cpu_seconds']['mean']:.2f}" in alltext
assert ('completed-subset averages' if e2e['incomplete_cases'] else 'All 140 reactions finish') in alltext.replace('\n',' ')
result=dict(status='passed',pages=len(pages),figures=figs,source_checks=True,references_resolved=True,no_overfull_boxes=True,
 manual_visual_review_required=True,scope='Numerical and build validation; visual review is recorded separately. MAPPA and the SLAP sweep use final-source campaigns; archived default-comparator outputs were rescored with the current evaluator. Completeness and limits are recorded in the evidence.',
 manuscript_sha256=hashlib.sha256((MAN/'manuscript.pdf').read_bytes()).hexdigest())
(MAN/'build/artifact-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
