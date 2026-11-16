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
tex=(MAN/'main.tex').read_text()+'\n'+'\n'.join(p.read_text() for p in (MAN/'includes').glob('*.tex') if p.name in ['paper.tex','supplement.tex','include-abstract.tex'])
for obsolete in ['Adaptive no-sweep','Separate experimental versions','Holdout follow-ups','earlier publication engine','N_{\\mathrm{order\\ changed}}']:
 assert obsolete not in tex,obsolete
assert not (MAN/'preprint.tex').exists()
assert not (MAN/'manuscript_bundle.zip').exists()
assert 'catalogue coverage' not in tex.lower()
for retired in ['generated-results', 'generated-cap-ablation-table', 'generated-coordinate-minima-table', 'Archived September']:
 assert retired not in tex, retired
for include in re.findall(r'\\input\{([^}]+)\}', tex):
 path=MAN/include
 assert path.with_suffix('.tex').is_file(), path
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
assert 'Archived' not in supplement and 'September' not in supplement
assert supplement.count(r'\section{') == 2
assert r'\label{tab:benchmark-settings}' in supplement
assert 'figs/fig2_golden.pdf' not in paper and 'figs/fig3_coordinate.pdf' not in paper
assert 'Final representation and complete window decoding' not in paper
bib=(MAN/'references.bib').read_text()
for group in re.findall(r'\\cite\w*\{([^}]+)\}',tex):
 for key in group.split(','):assert re.search(r'@\w+\{'+re.escape(key)+',',bib),key
figs=re.findall(r'\\includegraphics(?:\[[^]]*\])?\{(figs/[^}]+)\}',tex)
assert set(figs)=={'figs/fig1_algorithm.pdf','figs/fig4_alternatives.pdf','figs/fig6_current_benchmarks.pdf','figs/fig7_reference_offsets.pdf','figs/fig8_design_examples.pdf','figs/fig9_paired_slap.pdf','figs/fig10_minimum_choices.pdf','figs/fig11_saturation.pdf','figs/fig12_reference_percentages.pdf'}
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
for value in ['1,794','1,620','96.92',f"{choices['mean']:.2f}",f"{choices['single_choice']:,}"]:
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
assert 'first-returned' not in alltext.lower()
main_comparison=(MAN/'includes/generated-method-coverage-table.tex').read_text()
assert 'Unresolved' not in main_comparison
assert main_comparison.count('Single mapping')==5
for count in [1833,1842,1794,1620,1591,1547,1605,1561,692,1030]:
 assert f'{count:,} ({100*count/1851:.2f}' in main_comparison
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
result=dict(status='passed',pages=len(pages),figures=figs,source_checks=True,references_resolved=True,no_overfull_boxes=True,
 manual_visual_review_required=True,scope='Numerical and build validation; visual review is recorded separately. Checks cover the current Golden/EnzymeMap results, paired comparator timings and retained appendix only. The offset figure subtracts the saved detected minimum from a verified GT-consistent mapping count, including saved anchored completions. Every reaction contributes once, without additional optimization or lower-bound substitution.',
 manuscript_sha256=hashlib.sha256((MAN/'manuscript.pdf').read_bytes()).hexdigest())
(MAN/'build/artifact-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
