"""Current benchmark figures and tables from publication-safe saved evidence."""
import gzip
import json
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, MaxNLocator, NullFormatter


def evidence(man):
    summary = json.loads((man/'evidence/current_benchmarks.json').read_text())
    offsets = json.loads(gzip.decompress((man/'evidence/branch_offsets.json.gz').read_bytes()))
    for study in summary['studies']:
        assert study['directions'] == 2*study['reactions']
        for field in ('catalogue', 'decoded'):
            assert sum(study[field].values()) == study['reactions']
        assert study['config']['iso_tolerance'] == 2
        assert study['config']['branch_limit'] == 100
        assert study['decode_config']['event_slack'] == 5
        assert all(value > 0 for value in study['cpu_seconds'].values())
    rows = offsets['results']
    assert len(rows) == len({(r['dataset'], r['case']) for r in rows}) == 49825
    assert Counter(r['dataset'] for r in rows) == {'Golden': 1851, 'EnzymeMap': 47974}
    for row in rows:
        if row['status'] == 'complete':
            assert row['difference'] == row['branch_minimum']-row['overall_minimum']
        else:
            assert row.get('difference') is None
    return summary, offsets


def build(man):
    summary, offsets = evidence(man)
    studies = summary['studies']
    labels = [f"{s['dataset']} · {s['seeds']} seed"+('s' if s['seeds'] > 1 else '') for s in studies]
    colors = ['#087F8C', '#3366A5', '#77559B']

    def save(fig, stem):
        for ext in ('pdf', 'png', 'svg'):
            fig.savefig(man/'figs'/f'{stem}.{ext}', dpi=240, bbox_inches='tight', pad_inches=.08)
        plt.close(fig)

    fig, (coverage, cost) = plt.subplots(1, 2, figsize=(10.2, 4.1), layout='constrained',
                                        gridspec_kw={'width_ratios': [1.35, 1]})
    for y, study in enumerate(studies):
        found = study['decoded']['recovered']
        percent = 100*found/study['reactions']
        coverage.barh(y, percent, height=.36, color=colors[0])
        coverage.text(2, y, f"{found:,} / {study['reactions']:,}  ({percent:.3f}%)",
                      color='white', weight='bold', va='center', fontsize=9)
        coverage.text(1, y+.29, f"{study['decoded']['not_recovered']:,} absent; "
                      f"{study['decoded']['unknown']:,} unresolved", fontsize=8, color='#586974')
    coverage.set(yticks=range(3), yticklabels=labels, xlim=(0, 102), ylim=(2.65, -.5),
                 xlabel='Confirmed reactions (%)', xticks=[0, 25, 50, 75, 100])
    coverage.set_title('a  Decoded reference coverage', loc='left', weight='bold', pad=12)
    left = np.zeros(3)
    for stage, label, color in zip(('matching', 'deduplication', 'decoding'),
                                   ('Matching', 'Reduction', 'Decoding'), colors):
        values = np.array([s['cpu_seconds'][stage]/s['reactions'] for s in studies])
        cost.barh(range(3), values, left=left, height=.48, color=color, label=label)
        left += values
    for y, value in enumerate(left):
        cost.text(value+.25, y, f'{value:.2f}', va='center', fontsize=9)
    cost.set(yticks=range(3), yticklabels=labels, xlim=(0, max(left)*1.18), ylim=(2.65, -.8),
             xlabel='Mean compute CPU seconds / reaction')
    cost.set_title('b  Separate scientific stages', loc='left', weight='bold', pad=12)
    cost.legend(loc='upper left', ncol=3, fontsize=7, frameon=False)
    for ax in (coverage, cost):
        ax.grid(axis='x', alpha=.15)
        ax.set_axisbelow(True)
    save(fig, 'fig6_current_benchmarks')

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), layout='constrained')
    for ax, name, color in zip(axes, ('Golden', 'EnzymeMap'), colors[:2]):
        cohort = [r for r in offsets['results'] if r['dataset'] == name]
        counts = Counter(r['difference'] for r in cohort if r['status'] == 'complete')
        states = Counter(r['status'] for r in cohort)
        xs = sorted(counts)
        ax.bar(xs, [counts[x] for x in xs], width=.8, color=color)
        for x in xs:
            ax.text(x, counts[x]*1.15, f'{counts[x]:,}', ha='center', va='bottom', fontsize=8)
        ax.set_yscale('log')
        ax.set(xlim=(min(xs, default=0)-.7, max(5, max(xs, default=0))+.7),
               ylim=(.7, max(counts.values())*3), xlabel=r'$n_{\mathrm{GT}}$: branch minimum − detected minimum',
               ylabel='Reactions (log scale)')
        ax.set_title(f"{name}: {states['complete']:,} / {len(cohort):,} offsets", loc='left', weight='bold')
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
        ax.text(.5, -.23, f"Pending {states['pending']:,} · GT absent {states['not_recovered']:,} · "
                f"Unresolved {states['unknown']:,}", transform=ax.transAxes, ha='center', fontsize=8)
    fig.suptitle('Preliminary branch-event distribution · established offsets only', weight='bold', fontsize=11)
    save(fig, 'fig7_reference_offsets')

    coverage_rows, timing_rows = [], []
    for study in studies:
        name = f"{study['dataset']}, {study['seeds']} seed"+('s' if study['seeds'] > 1 else '')
        def rate(count):
            return f"{count:,} ({100*count/study['reactions']:.3f}\\%)"
        coverage_rows.append(' & '.join([name, f"{study['reactions']:,}",
            rate(study['catalogue']['recovered']), rate(study['decoded']['recovered']),
            str(study['decoded']['not_recovered']), str(study['decoded']['unknown'])])+r'\\')
        matching, reduction, decoding = [study['cpu_seconds'][s]/study['reactions']
                                         for s in ('matching', 'deduplication', 'decoding')]
        timing_rows.append(' & '.join([name, *[f'{v:.2f}' for v in
            (matching, reduction, matching+reduction, decoding, matching+reduction+decoding)]])+r'\\')
    def table(name, header, rows):
        text = r'\begin{tabular}{lrrrrr}\toprule'+'\n'+header+r'\\\midrule'+'\n'
        text += '\n'.join(rows)+'\n'+r'\bottomrule\end{tabular}'+'\n'
        (man/'includes'/name).write_text(text)
    table('generated-current-coverage-table.tex',
          'Study & Reactions & Catalogue & Decoded & Absent & Unknown', coverage_rows)
    table('generated-current-timing-table.tex',
          'Study & Matching & Reduction & Mapping & Decoding & Total', timing_rows)
