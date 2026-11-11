"""Validate the final manuscript's numerical sources and compiled references."""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1]/"build/mpl-cache"))
import gzip,hashlib,json,re
from publication_privacy import check_publication
from pypdf import PdfReader
MAN=Path(__file__).resolve().parents[1];E=MAN/'evidence'
check_publication('manuscript.pdf',(MAN/'manuscript.pdf').read_bytes())
def read(n):return json.loads((E/n).read_text())
# Separate current catalogue and decoded coverage; never turn unknowns into misses.
from build_benchmark_figures import anchored_evidence, evidence
current, offsets = evidence(MAN)
studies = {row['id']: row for row in current['studies']}
assert [studies[k]['decoded']['recovered'] for k in ('golden1','golden10','enzyme1')] == [1833,1842,47796]
assert [studies[k]['catalogue']['recovered'] for k in ('golden1','golden10','enzyme1')] == [1834,1843,47798]
assert [studies[k]['decoded']['unknown'] for k in ('golden1','golden10','enzyme1')] == [0,0,0]
closure=E
final=json.loads((closure/'decoded_verification_closure.json').read_text())
directions=json.loads(gzip.decompress((closure/'decoded_verification_closure_directions.json.gz').read_bytes()))
assert final['source_hashes_unchanged'] and not final['core_algorithm_changed']
assert all(proof['valid'] for proof in final['independent_replays'])
assert hashlib.sha256((closure/'decoded_verification_closure_directions.json.gz').read_bytes()).hexdigest()==final['direction_snapshot_sha256']
for arm,study in [('seeds1','golden1'),('seeds10','golden10')]:
 rows=directions['studies'][arm]['reactions']
 assert len(rows)==1851 and len({r['case'] for r in rows})==1851
 assert final['studies'][arm]['reaction_counts']==studies[study]['decoded']
 for row in rows:
  states=row['directions'].values()
  assert row['status']==('recovered' if 'recovered' in states else 'not_recovered')
  if row['status']=='not_recovered':assert set(states)=={'not_recovered'}
 assert sum(r['status']=='recovered' for r in rows)==studies[study]['decoded']['recovered']
assert offsets['additional_searches'] == 0
from build_benchmark_figures import offset_counts
for dataset, study in [('Golden', 'golden1'), ('EnzymeMap', 'enzyme1')]:
 rows=[r for r in offsets['results'] if r['dataset']==dataset]
 assert sum(r['source']=='verified' for r in rows)==studies[study]['decoded']['recovered']
 assert sum(offset_counts(rows).values())==studies[study]['reactions']
from build_benchmark_figures import comparison_evidence
paired, choices, choice_rows = comparison_evidence(MAN)
anchored = anchored_evidence(MAN)
assert hashlib.sha256((E/'branch_offsets.json.gz').read_bytes()).hexdigest() == anchored['original_branch_offsets_sha256']
for key, filename in [('baseline','branch_offsets.json.gz'),
                      ('golden','decoded_verification_closure_directions.json.gz'),
                      ('anchored','anchored_offsets.json')]:
 assert hashlib.sha256((E/filename).read_bytes()).hexdigest()==offsets['source_sha256'][key]
from collections import Counter
expected_anchored_counts = {
 'Golden': {0:9, 2:2, 3:1, 4:2, 5:2, 7:1, 12:1},
 'EnzymeMap': {-2:9, 0:83, 2:32, 4:51, 6:3},
}
for dataset, expected in expected_anchored_counts.items():
 assert Counter(r['difference'] for r in anchored['results'] if r['dataset']==dataset) == expected
expected_merged_counts = {
 'Golden': {0:1630, 1:90, 2:88, 3:24, 4:11, 5:6, 7:1, 12:1},
 'EnzymeMap': {-2:9, 0:44086, 1:1384, 2:2190, 3:122, 4:171, 5:9, 6:3},
}
for dataset, expected in expected_merged_counts.items():
 rows=[r for r in offsets['results'] if r['dataset']==dataset]
 merged=offset_counts(rows)
 assert merged==expected and sum(merged.values())==len(rows)
original={(r['dataset'],r['case']):r for r in json.loads(gzip.decompress((E/'branch_offsets.json.gz').read_bytes()))['results']}
completed={(r['dataset'],r['case']):r for r in anchored['results']}
certificates={}
for row in directions['studies']['seeds1']['directions']:
 if row['status']=='recovered':certificates.setdefault(row['case'],[]).append(row)
for row in offsets['results']:
 key=row['dataset'],row['case']
 if row['source']=='anchored':
  source=completed[key]
  assert original[key]['status']=='not_recovered'
  assert row['reference_count']==source['anchored_count']
 else:
  source=original[key]
  assert source['status']!='not_recovered'
  if row['dataset']=='Golden':
   certificate=min(certificates[row['case']],key=lambda r:(r['event_count'],r['reversed']))
   assert row['reference_count']==certificate['event_count']
   assert row['certificate_sha256']==certificate['certificate_sha256']
 assert row['overall_minimum']==source['overall_minimum']
example=next(r for r in offsets['results'] if r['dataset']=='Golden' and r['case']==268)
assert (example['reference_count'],example['overall_minimum'],example['difference'])==(16,12,4)
golden_minima = {r['case']:r['detected_minimum'] for r in choice_rows}
assert all(r['overall_minimum']==golden_minima[r['case']] for r in anchored['results'] if r['dataset']=='Golden')
# Reference-supplied completions cannot become choices returned without anchors.
historical_choices = json.loads(gzip.decompress((E/'minimum_choices_reactions.json.gz').read_bytes()))
choice_membership = {r['case']:r['reference_at_minimum'] for r in historical_choices}
assert all(choice_membership[r['case']] is False for r in anchored['results'] if r['dataset']=='Golden')
assert sum(r['reference_reindexed'] for r in anchored['results']) == 11
assert paired['coverage']['mappa']['recovered'] == 1839
assert paired['coverage']['slap_sweep']['recovered'] == 1794
assert paired['coverage']['slap_bidirectional']['recovered'] == 1620
assert choices['additional_searches'] == 0
original_branches = {r['case']: r for r in json.loads(gzip.decompress((E/'minimum_branches_reactions.json.gz').read_bytes()))}
for row in choice_rows:
 original = original_branches[row['case']]
 assert all(row[key] == original[key] for key in ('detected_minimum', 'branches'))
 assert sorted(row['decoded_hashes']) == sorted(original['decoded_hashes'])
 assert 1 <= row['reaction_cores'] <= row['saved_mappings']
 assert row['reaction_cores'] == len({center['id'] for center in row['centers']})
 assert row['reaction_cores'] == sum(row['physical_event_counts'].values())
from build_molecule_figure import check_example, check_design_examples
assert read('molecule_example.json') == json.loads(json.dumps(check_example()))
assert read('design_examples.json') == check_design_examples()
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
assert sources['historical_campaign']['fresh_campaign_complete']
for row in sources['snapshots']:
 assert hashlib.sha256((E/row['snapshot']).read_bytes()).hexdigest()==row['sha256'],row['snapshot']
tex=(MAN/'main.tex').read_text()+'\n'+'\n'.join(p.read_text() for p in (MAN/'includes').glob('*.tex') if p.name in ['paper.tex','supplement.tex','include-abstract.tex'])
for obsolete in ['Adaptive no-sweep','Separate experimental versions','Holdout follow-ups','earlier publication engine','N_{\\mathrm{order\\ changed}}']:
 assert obsolete not in tex,obsolete
assert not (MAN/'preprint.tex').exists()
assert not (MAN/'manuscript_bundle.zip').exists()
assert 'catalogue coverage' not in tex.lower()
coverage_table=(MAN/'includes/generated-current-coverage-table.tex').read_text()
assert 'Catalogue' not in coverage_table
for study in current['studies']:
 assert f"{study['decoded']['recovered']:,}" in coverage_table
paper=(MAN/'includes/paper.tex').read_text()
supplement=(MAN/'includes/supplement.tex').read_text()
stages=['matching-rules','growth','completion','branching','cut-sweep','reduction','decoding','verification']
positions=[paper.index(r'\label{sec:'+stage+'}') for stage in stages]
assert positions==sorted(positions)
assert r'\subsubsection{Growth matching:' in paper
assert r'\subsubsection{Saturation matching:' in paper
assert paper.index(r'\label{sec:completion}') < paper.index(r'\label{fig:saturation}') < paper.index(r'\label{sec:branching}')
assert r'\subsubsection{Hierarchical branching:' in paper
assert r'\label{tab:growth-transitions}' in paper
for transition in [r'0\to0', r'1\to0', r'\mathrm{many}\to0', r'1\to1',
                   r'1\to\mathrm{many}', r'\mathrm{many}\to1', r'\mathrm{many}\to\mathrm{many}']:
 assert transition in paper,transition
assert 'One failed edge does not end growth' in paper
assert paper.index(r'\label{sec:verification}') > paper.index(r'\section{Evaluation}')
assert 'fig5_pipeline' not in tex and 'fig:pipeline' not in tex
assert paper.index(r'\label{fig:alternatives}') < paper.index(r'\section{Methods}')
assert 'Reduction Stage' not in tex and 'Three stages' not in paper
assert paper.index(r'\label{tab:method-coverage}') > paper.index(r'\section{Results}')
assert 'singleton' in paper
assert not re.search(r'\b(?:reactants?|products?)\b',tex,re.I)
assert r'\emph{base}' in paper and r'\emph{target}' in paper
assert 'V_B' in supplement and 'V_T' in supplement
assert not re.search(r'\b[VZW]_[RP]\b',tex)
assert all(study['catalogue']['unknown'] == 0 for study in current['studies'])
assert r'\label{sec:historical-evaluation}' in supplement
assert 'figs/fig2_golden.pdf' not in paper and 'figs/fig3_coordinate.pdf' not in paper
assert 'Final representation and complete window decoding' not in paper
bib=(MAN/'references.bib').read_text()
for group in re.findall(r'\\cite\w*\{([^}]+)\}',tex):
 for key in group.split(','):assert re.search(r'@\w+\{'+re.escape(key)+',',bib),key
figs=re.findall(r'\\includegraphics(?:\[[^]]*\])?\{(figs/[^}]+)\}',tex)
assert set(figs)=={'figs/fig1_algorithm.pdf','figs/fig2_golden.pdf','figs/fig3_coordinate.pdf','figs/fig4_alternatives.pdf','figs/fig6_current_benchmarks.pdf','figs/fig7_reference_offsets.pdf','figs/fig8_design_examples.pdf','figs/fig9_paired_slap.pdf','figs/fig10_minimum_choices.pdf','figs/fig11_saturation.pdf','figs/fig12_reference_percentages.pdf'}
assert 'figs/fig12_reference_percentages.pdf' in supplement
assert 'figs/fig12_reference_percentages.pdf' not in paper
for fig in figs:
 p=MAN/fig
 assert len(PdfReader(p).pages)==1
 for ext in ['svg','png']:assert p.with_suffix('.'+ext).is_file()
log=(MAN/'build/main.log').read_text()
assert not re.search(r'(?:Citation|Reference).*undefined|There were undefined|Overfull \\[hv]box|Missing character:',log), 'Inspect TeX log'
r=PdfReader(MAN/'manuscript.pdf');pages=[p.extract_text() for p in r.pages];alltext='\n'.join(pages)
assert '??' not in alltext
for value in ['1,833','1,842','47,796','99.629','49,825','1,630','44,086','88.06','91.90']:
 assert value in alltext,value
for value in ['1,794','1,620','96.921',f"{choices['mean']:.2f}",f"{choices['single_choice']:,}"]:
 assert value in alltext,value
assert 'cache control' not in paper.lower() and 'cache speedup' not in paper.lower()
assert 'approximately linearly' in paper
assert 'binary and weighted' in paper and 'single-edge cut sweep' in paper
for study in current['studies']:
 means=[study['cpu_seconds'][s]/study['reactions'] for s in ('matching','deduplication','decoding')]
 for value in [sum(means[:2]),means[2],sum(means)]:
  assert f'{value:.2f}' in alltext
assert 'ordinary-boundary' in tex and 'starting orbit' in tex
results=paper.split(r'\section{Results}',1)[1].split(r'\section{Discussion}',1)[0]
assert results.count(r'\subsection{')==4
assert not any(word in results.lower() for word in ['pending', 'preliminary', 'first-returned'])
assert 'generated-method-coverage-table' in results and 'ten seeds' in results.lower()
assert 'lower endpoint' not in results and 'no additional optimization' in results
assert 'Saved anchored completions' in results
assert 'capped lines' not in results and 'upper count bounds' not in results
assert 'whole-growth replay' not in paper
assert not re.search(r'\bMac\b|same-host|worker count|ran concurrently',alltext), 'Keep execution details in the reproducibility records'
for value in [f"{seed['methods']['seeds1']['golden_recovery_percent']:.2f}",f"{seed['methods']['seeds3']['golden_recovery_percent']:.2f}",'120,052','300','162','168']:
 assert value in alltext,value
assert 'Golden reference recovery and runtime' in alltext
assert 'First correct' not in alltext and 'two ten-order cases remain unresolved' not in alltext
assert 'first-returned' not in alltext.lower()
main_comparison=(MAN/'includes/generated-method-coverage-table.tex').read_text()
assert 'Unresolved' not in main_comparison
assert main_comparison.count('Single mapping')==5
for count in [1833,1842,1794,1620,1591,1547,1605,1561,692,1030]:
 assert f'{count:,} ({100*count/1851:.2f}' in main_comparison
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
 figure_text=PdfReader(MAN/fig).pages[0].extract_text()
 if fig.endswith('fig6_current_benchmarks.pdf'):
  for label in ['99.028','99.514','99.629','Matching','Decoding']:
   assert label in figure_text,label
 elif fig.endswith('fig7_reference_offsets.pdf'):
  for label in ['Golden','EnzymeMap','1,851','47,974',
                '1,630','44,086','Event offset n']:
   assert label in figure_text,label
  assert not any(label in figure_text for label in ['Pending','GT absent','Possible upper count','Anchored completion','Retained branches'])
 elif fig.endswith('fig12_reference_percentages.pdf'):
  for label in ['Golden','EnzymeMap','1,851','47,974','88.06%','91.90%','Reactions (%)']:
   assert label in figure_text,label
  for dataset in ('Golden','EnzymeMap'):
   cohort=[r for r in offsets['results'] if r['dataset']==dataset]
   percentages=[100*count/len(cohort) for count in offset_counts(cohort).values()]
   assert abs(sum(percentages)-100)<1e-10
   for value in percentages:
    assert (f'{value:.2f}%' if value>=.1 else f'{value:.3f}%') in figure_text
 elif fig.endswith('fig8_design_examples.pdf'):
  for label in ['Greedy constraints','1 inner + 0 outer','1 inner + 2 outer','2 outer events','Base','Target']:
   assert label in figure_text,label
 else:
  assert 'MAPPA' in figure_text,fig
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
assert not re.search(r'\b(?:competition mode|takeover decoding)\b',tex,re.I)
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
assert 'Archived MAPPA search plus event-window queries' in alltext.replace('\n',' ')
assert f"{e2e['complete_cohort']['end_to_end']['cpu_seconds']['mean']:.2f}" in alltext
assert ('completed-subset averages' if e2e['incomplete_cases'] else 'All 140 reactions finish') in alltext.replace('\n',' ')
result=dict(status='passed',pages=len(pages),figures=figs,source_checks=True,references_resolved=True,no_overfull_boxes=True,
 manual_visual_review_required=True,scope='Numerical and build validation; visual review is recorded separately. Current Golden and EnzymeMap coverage and timings are validated separately from archived comparisons. The offset figure subtracts the saved detected minimum from a verified GT-consistent mapping count, including saved anchored completions. Every reaction contributes once, without additional optimization or lower-bound substitution.',
 manuscript_sha256=hashlib.sha256((MAN/'manuscript.pdf').read_bytes()).hexdigest())
(MAN/'build/artifact-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
