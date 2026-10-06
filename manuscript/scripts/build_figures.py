"""Build the final paper's vector figures from checked-in evidence; no AAM runs."""
import os,json
from pathlib import Path
MAN=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(MAN/'build/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':240})
INK='#223340';MUTED='#586974';GREEN='#087F8C';BLUE='#3366A5';ORANGE='#CB702A';RED='#B63F52';LIGHT='#D5DEE2';PURPLE='#77559B'
E=MAN/'evidence';F=MAN/'figs';F.mkdir(exist_ok=True)
def read(n):return json.loads((E/n).read_text())
def save(fig,n):
 for ext in ['pdf','svg','png']:fig.savefig(F/f'{n}.{ext}',facecolor='white',bbox_inches='tight',pad_inches=.05)
 plt.close(fig)
from build_molecule_figure import build as build_molecule_figure
build_molecule_figure(MAN)
from build_alternatives_figure import build as build_alternatives_figure
build_alternatives_figure(MAN)

timing=read('timing_comparison.json')
competitors=read('competitors.json')
local=next(d for d in competitors['methods'] if d['method']=='localmapper')
seed=read('seed_comparison.json');methods=seed['methods'];slap=read('slap_sweep.json');N=1851
keys=['seeds1','seeds2','seeds3','seeds10'];counts=[methods[k]['golden_outcomes']['recovered'] for k in keys];costs=[methods[k]['common_mean_cpu_seconds'] for k in keys]
assert seed['fresh'] and all(sum(methods[k]['golden_outcomes'].values()) == N for k in keys)
assert len(seed['common_case_indices']) > 0
fig,axs=plt.subplots(1,2,figsize=(9.2,4.7),gridspec_kw={'width_ratios':[1.4,1.3]},layout='constrained')
a=axs[0];labels=['MAPPA · 1 seed\n(default)','MAPPA · 2 seeds','MAPPA · 3 seeds','MAPPA · 10 seeds','SLAP sweep\n(prior baseline)','LocalMapper\n(prior SOTA, 2024)'];vals=[*counts,slap['sweep_union_recovered'],local['any_correct']];y=np.arange(6)
a.barh(y,np.array(vals)/N*100,color=[GREEN,BLUE,'#52789F',PURPLE,ORANGE,MUTED],height=.58)
for i,v in enumerate(vals):a.text(v/N*100-1.5,i,f'{v:,} / {N:,}  ({v/N*100:.2f}%)',ha='right',va='center',color='white',fontsize=9,weight='bold')
a.set(yticks=y,yticklabels=labels,xlim=(0,102),xticks=[0,25,50,75,100],xlabel='Reference recovery over returned output (%)');a.invert_yaxis();a.set_title('a  MAPPA search-stage recovery',loc='left',weight='bold',pad=13)
a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
shown=[costs[0],costs[1],seed['same_host_timing']['metrics']['slap']['mean_seconds']]
a=axs[1];ys=np.arange(3);policies=['smaller_first','larger_first','bidirectional']
for j,(prefix,label,color) in enumerate([('mappa1','MAPPA · 1 seed (default)',GREEN),('mappa2','MAPPA · 2 seeds',BLUE),('slap','SLAP sweep (prior baseline)',ORANGE)]):
 vals=[timing['methods'][prefix+'_'+p]['stats']['mean'] for p in policies];positions=ys+(j-1)*.23
 a.barh(positions,vals,color=color,height=.21,label=label)
 for y,v in zip(positions,vals):a.text(v+.06,y,f'{v:.2f}',va='center',fontsize=8)
a.set(yticks=ys,yticklabels=['Smaller-first','Larger-first','Bidirectional'],xlim=(0,6.2),ylim=(2.65,-1.2),xlabel='Mean CPU seconds / reaction')
a.set_title('b  Search cost by size orientation',loc='left',weight='bold',pad=13)
a.legend(loc='upper left',frameon=False,fontsize=8);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
save(fig,'fig2_golden')
uncut=read('unswept.json')

dedup=read('final_dedup.json');rows=dedup['per_case']
assert len(rows)==140 and sum(r['flat_saved_families'] for r in rows)==dedup['flat_families']
assert dedup['competition'] is False and all(r['complete'] for r in rows)
coord=read('coordinate_minima.json')
fig,axs=plt.subplots(2,2,figsize=(8.8,6.1),layout='constrained',gridspec_kw={'hspace':.16,'wspace':.10})
cs=[coord['comparisons'][k] for k in ['slap_sweep','native_slap']]
a=axs[0,0]
for i,d in enumerate(cs):
 low=d['count_relations'].get('lower',0);eq=d['count_relations'].get('equal',0);assert low+eq==140
 a.barh(i,low,color=GREEN,height=.45,label='MAPPA fewer' if i==0 else None)
 a.barh(i,eq,left=low,color=BLUE,height=.45,label='Equal count' if i==0 else None)
 a.text(low/2,i,str(low),ha='center',va='center',color='white',weight='bold',fontsize=8)
 a.text(low+eq/2,i,str(eq),ha='center',va='center',color='white',weight='bold')
a.set(yticks=[0,1],yticklabels=['vs SLAP sweep','vs native SLAP'],xlim=(0,145),ylim=(1.85,-.8),xlabel='Reactions (140 total)')
a.set_title('a  MAPPA minimum event counts',loc='left',weight='bold',pad=12);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a.legend(loc='upper left',fontsize=8,frameon=False,ncol=2);a.text(.02,.06,'MAPPA has a higher minimum in 0 reactions.',transform=a.transAxes,fontsize=8,color=MUTED)
a=axs[0,1]
for i,d in enumerate(cs):
 n=d['shared_minimum_patterns'];total=d['equal_minimum_patterns'];a.barh(i,n,color=GREEN,height=.45);a.barh(i,total-n,left=n,color=RED,height=.45)
 a.text(n/2,i,f'{n}/{total}',ha='center',va='center',color='white',weight='bold')
a.set(yticks=[0,1],yticklabels=['SLAP sweep','Native SLAP'],xlim=(0,175),ylim=(1.85,-.6),xlabel='Comparator patterns retained by MAPPA')
a.set_title('b  Patterns where minima agree',loc='left',weight='bold',pad=12);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a.text(.02,.03,'136 equal-minimum cases for the sweep;\n125 for native SLAP. Red: missing alternatives.',transform=a.transAxes,fontsize=7.8,color=MUTED)
a=axs[1,0]
series=[('MAPPA',GREEN,[r['mappa_minimum'] for r in coord['per_case']]),
        ('SLAP sweep',ORANGE,[r['comparators']['slap_sweep']['minimum'] for r in coord['per_case']]),
        ('Native SLAP',MUTED,[r['comparators']['native_slap']['minimum'] for r in coord['per_case']])]
values=list(range(max(max(v) for _,_,v in series)+1));largest=0
for j,(label,color,mins) in enumerate(series):
 freq=[mins.count(k) for k in values];largest=max(largest,max(freq))
 a.bar(np.array(values)+(j-1)*.25,freq,color=color,width=.23,label=label)
a.set(xticks=values,ylim=(0,largest*1.35),xlabel='Minimum signed-event count',ylabel='Reactions')
a.set_title('c  Minimum-event distributions',loc='left',weight='bold',pad=12);a.grid(axis='y',alpha=.15);a.set_axisbelow(True)
a.legend(loc='upper left',ncol=3,fontsize=7,frameon=False)
a=axs[1,1];keys_stage=[k for k in ['mappa_search','mappa_decoding','slap_mapping','slap_h_refinement'] if k in coord['timing']];means=[coord['timing'][k]['stats']['mean'] for k in keys_stage];ys=np.arange(len(keys_stage))
a.barh(ys,means,color=[GREEN if k=='mappa_search' else PURPLE if k=='mappa_decoding' else ORANGE for k in keys_stage],height=.52)
for i,v in enumerate(means):a.text(v+.12,i,f'{v:.2f}',va='center',fontsize=8)
a.set(yticks=ys,yticklabels=[dict(mappa_search='MAPPA search',mappa_decoding='MAPPA window decode',slap_mapping='SLAP sweep mapping',slap_h_refinement='SLAP H refinement')[k] for k in keys_stage],xlim=(0,max(means)*1.25),ylim=(len(keys_stage)-.3,-.7),xlabel='Mean CPU seconds / reaction')
a.set_title('d  Recorded stage costs',loc='left',weight='bold',pad=12);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
save(fig,'fig3_coordinate')
windows=sorted((int(k),v) for k,v in dedup['window_distribution'].items())
(MAN/'includes/generated-decoder-table.tex').write_text('\\begin{tabular}{l'+ 'r'*len(windows)+'}\\toprule\nMaximum events & '+' & '.join(str(k) for k,v in windows)+r'\\'+'\nReactions (all complete) & '+' & '.join(str(v) for k,v in windows)+r'\\'+'\n'+r'\bottomrule\end{tabular}'+'\n')
paired=slap['paired_recovery']['seeds1']
macros=dict(GoldenOneRecovered=f"{counts[0]:,}",GoldenOnePercent=f"{100*counts[0]/N:.2f}",
 GoldenThreeRecovered=f"{counts[2]:,}",GoldenThreePercent=f"{100*counts[2]/N:.2f}",
 GoldenTenRecovered=f"{counts[3]:,}",GoldenTenPercent=f"{100*counts[3]/N:.2f}",GoldenTenUnknown=str(methods["seeds10"]["golden_outcomes"]["unknown"]),SlapRecovered=f"{slap['sweep_union_recovered']:,}",
 SlapPercent=f"{100*slap['sweep_union_recovered']/N:.2f}",
 RecoveryDifference=f"{100*(counts[0]-slap['sweep_union_recovered'])/N:.2f}",
 CommonCases=f"{len(seed['common_case_indices']):,}",OneCPU=f"{costs[0]:.2f}",TwoCPU=f"{costs[1]:.2f}",SlapCPU=f"{shown[2]:.2f}",
 AamOnly=str(len(paired['aam_only'])),SlapOnly=str(len(paired['slap_only'])),
 SweepRecovered=str(coord['comparisons']['slap_sweep']['covered_in_mappa_window']),
 SweepClasses=str(coord['comparisons']['slap_sweep']['minimum_patterns']),
 SlapAbsent=str(slap['outcomes']['not_recovered']),SlapUnknown=str(slap['outcomes']['unknown']),
 FlatFamilies=f"{dedup['flat_families']:,}",DecodedClasses=str(dedup['event_classes']))
(MAN/'includes/generated-results.tex').write_text('% Generated from complete fresh benchmark evidence.\n'+''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in macros.items()))
print('Built four main figures and complete-campaign tables and numerical macros.')

# One main comparison: primary recovery plus explicitly identified timing cohorts.
competitors=read('competitors.json')
assert competitors['rescored_with_current_evaluator'] and competitors['unchanged_case_outcomes']
controlled=read('golden_cap_ablation.json')
controlled_rows={r['key']:r for r in controlled['rows']}
call_times={r['key']:r for r in timing['default_comparators']}
labels={'rxnmapper':'RXNMapper','localmapper':r'LocalMapper$^{\dagger}$','chython':'Chython',
        'slap_binary':'SLAP, binary','slap_weighted':'SLAP, weighted','indigo':'Indigo','rdt':'RDT'}
rows=[]
audit=[]
def comparison_row(name,setting,output,n,cpu,key,timing_group):
 assert 0 <= n <= N and cpu >= 0
 audit.append(dict(method=name,setting=setting,output=output,recovered=n,denominator=N,
                   mean_cpu_seconds=cpu,timing_key=key,timing_group=timing_group))
 marker=r"$^{a}$" if timing_group=='B' else ''
 return f"{name} & {setting} & {output} & {n:,} ({100*n/N:.2f}\\%) & {cpu:.3f}{marker} " + r"\\"
for k,nseed in zip(keys,[1,2,3,10]):
 tkey=f'mappa_c100_s{nseed}_sweep';cost=controlled_rows[tkey]
 n=methods[k]['golden_outcomes']['recovered']
 assert len(cost['cases']['recovered'])==n
 rows.append(comparison_row('MAPPA',f'{nseed} seed'+('s' if nseed>1 else '')+', sweep'+(' (default)' if nseed==1 else ''),'Families',n,cost['paired_mean'],tkey,'A'))
tkey='mappa_c100_s1_uncut'
rows.append(comparison_row('MAPPA','1 seed, no sweep','Families',uncut['methods']['mappa']['counts']['recovered'],controlled_rows[tkey]['paired_mean'],tkey,'A'))
for key,setting,n in [('slap_uncut','No sweep',uncut['methods']['slap']['counts']['recovered']),('slap_sweep','Sweep',slap['sweep_union_recovered'])]:
 name=r'SLAP, bidirectional$^{\ddagger}$' if key=='slap_sweep' else 'SLAP, bidirectional'
 rows.append(comparison_row(name,setting,'Candidates',n,controlled_rows[key]['paired_mean'],key,'A'))
for key in ['slap_binary','rxnmapper','localmapper','chython','indigo','rdt']:
 d=next(d for d in competitors['methods'] if d['method']==key)
 assert d['total']==N and len(d['any_correct_cases'])==d['any_correct']
 rows.append(comparison_row(labels[key],'Default, no sweep','Candidates' if key.startswith('slap') else 'One bijection',d['any_correct'],call_times[key]['mean_cpu_seconds'],key,'B'))
(MAN/'includes/generated-competitor-table.tex').write_text(r"\begin{tabular}{@{}lllr r@{}}\toprule"+'\n'+r"Method & Search setting & Output & Reference recovered & Mean CPU s\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')
(MAN/'build/golden-comparison-table.json').write_text(json.dumps(dict(rows=audit,
 recovery_scope='Primary evaluation: all 1,851 records; failures and unknowns retained.',
 timing_A='golden_cap_ablation.json: paired_mean, common 1,403 completed reactions, AMD EPYC 9J14, search only.',
 timing_B='timing_comparison.json: default_comparators mean_cpu_seconds, completed calls including invalid outputs, CPU model unrecorded.',
 no_cross_group_speed_ranking=True),indent=2)+'\n')

# The case-level failure table is generated from the audited final misses.
misses=read('golden_miss_analysis.json')
labels={'extra_matched_pair':'One extra matched pair',
 'atom_correspondence_excluded_by_relaxed_orbits':'Atom pairing excluded by orbit bounds',
 'joint_correspondence_excluded_by_full_verifier':'Joint C/O correspondence excluded'}
rows=[f"{r['case']} & {r['reference_pairs']} & {r['retained_pair_counts'][0]} & {labels[r['classification']]} " + r"\\" for r in misses['rows']]
(MAN/'includes/generated-miss-table.tex').write_text(r"\begin{tabular}{lrrl}\toprule"+'\n'+r"Case & Reference pairs & Retained pairs & Observed mismatch\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')

# Pair orientation-specific coverage with matched-cohort runtime distributions.
rows=[]
for prefix,label in [('mappa1','MAPPA, 1 seed (default)'),('mappa2','MAPPA, 2 seeds'),('slap','SLAP sweep (prior baseline)')]:
 if rows:rows.append(r"\midrule")
 for i,(policy,dlabel) in enumerate([('smaller_first','Smaller-first'),('larger_first','Larger-first'),('bidirectional','Bidirectional')]):
  d=timing['methods'][prefix+'_'+policy];n=d['recovered'];st=d['stats']
  row=[label if i==0 else '',dlabel,f'{n:,} ({100*n/N:.2f}\\%)']+[f'{st[k]:.3f}' for k in ['mean','median','p95']]
  row.append(f"{st['mean']/timing['methods']['slap_'+policy]['stats']['mean']:.2f}")
  rows.append(' & '.join(row)+r"\\")
(MAN/'includes/generated-direction-table.tex').write_text(r"\begin{tabular}{@{}llrrrrr@{}}\toprule"+'\n'+r"Method & Orientation & Recovered & Mean & Median & 95th pct. & CPU/SLAP\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')
# Minimum-event comparison is separate from mapping-accuracy evaluation.
rows=[]
for name,label in [('slap_sweep','SLAP sweep'),('native_slap','Native SLAP')]:
 d=coord['comparisons'][name];n=d['shared_minimum_patterns'];total=d['equal_minimum_patterns']
 cells=[label]+[str(d['count_relations'].get(k,0)) for k in ['lower','equal','higher']]+[f'{n}/{total}']
 rows.append(' & '.join(cells)+r"\\")
(MAN/'includes/generated-coordinate-minima-table.tex').write_text(r"\begin{tabular}{@{}lrrrr@{}}\toprule"+'\n'+r"Compared with & MAPPA fewer & Equal & MAPPA more & Patterns at equal minima\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')
rows=[]
for key,label in [('mappa_search','MAPPA sweep search'),('mappa_decoding','MAPPA full-window decoding'),('slap_mapping','SLAP sweep mapping'),('slap_h_refinement','SLAP H refinement / initial scoring')]:
 if key not in coord['timing']:continue
 d=coord['timing'][key]['stats'];cells=[label]+[f'{d[k]:.3f}' for k in ['mean','median','p95']]
 rows.append(' & '.join(cells)+r"\\")
(MAN/'includes/generated-coordinate-time-table.tex').write_text(r"\begin{tabular}{@{}lrrr@{}}\toprule"+'\n'+r"Recorded stage & Mean CPU s & Median CPU s & 95th pct. CPU s\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')


# Final public pipeline: keep fresh end-to-end timings separate from archived stages.
e2e=read('end_to_end_timing.json')
assert e2e['attempted']==e2e['requested']==140 and e2e['all_completed_match']
n=e2e['completed'];missing=len(e2e['incomplete_cases']);st=e2e['complete_cohort']
rows=[]
for key,label in [('search','MAPPA sweep search'),('catalogue','Final-family catalogue'),('decode','Bond-event decoding'),('serialize','Candidate serialization'),('end_to_end','End to end')]:
 d=st[key]['cpu_seconds'];rows.append(' & '.join([label]+[f"{d[k]:.3f}" for k in ['mean','median','p95']])+r"\\")
(MAN/'includes/generated-end-to-end-table.tex').write_text(r"\begin{tabular}{@{}lrrr@{}}\toprule"+'\n'+r"Stage & Mean CPU s & Median CPU s & 95th pct. CPU s\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')
total=st['end_to_end']['cpu_seconds'];decode=st['decode']['cpu_seconds']
cohort = "completed subset" if missing else "full set"
text=(f"A fresh run of the final implementation completed search and event decoding for {n}/140 coordinate reactions. "
 f"On this {cohort}, the full pipeline takes a mean {total['mean']:.2f} CPU seconds per reaction "
 f"(median {total['median']:.2f}; 95th percentile {total['p95']:.2f}). Decoding accounts for a mean {decode['mean']:.2f} CPU seconds "
 f"(median {decode['median']:.2f}). "
 + (f"The remaining {missing} reactions reached the five-minute watchdog; completed-subset averages do not estimate full-dataset completion cost. " if missing else "")
 + r"These measurements include search, final-family construction, decoding, and candidate output; they are separate from the search-only comparison above (\Cref{tab:end-to-end})."+'\n')
(MAN/'includes/generated-end-to-end-text.tex').write_text(text)
wall=st['end_to_end']['wall_seconds'];cpu=e2e['manifest']['cpu']
scope=(f"We measured one-seed, cap-2,000 cut-sweep search followed by the public event decoder on the {cpu} CPU. "
 "CPU times are reported per reaction. "
 "The event thresholds and windows match the bond-event comparison, with no class-count cap. "
 f"All {n} completed reactions reproduce the archived event-class sets and family counts. "
 f"Median end-to-end wall time is {wall['median']:.2f} seconds. "
 "Timing includes input loading, checkpoint I/O, final-family construction, decoding, and candidate serialization; "
 "initial imports and reference-set verification are excluded. "
 + ("Interrupted cases are indexed "+', '.join(map(str,e2e['incomplete_cases']))+" in the released records; their partial results are not completion certificates. " if missing else "All 140 reactions finish within the five-minute watchdog. ")
 + "Detailed timing records and reproduction instructions are provided with the code.\n")
(MAN/'includes/generated-end-to-end-scope.tex').write_text(scope)
