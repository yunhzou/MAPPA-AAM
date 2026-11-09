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
    build_comparisons(man, save)


def comparison_evidence(man):
    paired = json.loads((man/'evidence/paired_slap.json').read_text())
    choices = json.loads((man/'evidence/minimum_choices.json').read_text())
    rows = json.loads(gzip.decompress((man/'evidence/minimum_choices_reactions.json.gz').read_bytes()))
    assert paired['reactions'] == choices['expected'] == len(rows) == 1851
    assert len({row['case'] for row in rows}) == 1851
    for result in paired['coverage'].values():
        assert sum(result.values()) == 1851
    assert paired['verification_audit']['disagreements'] == 0
    exact = [row for row in rows if row['status'] == 'complete']
    assert len(exact) == choices['complete_counts'] == 1829
    assert Counter(row['choices'] for row in exact) == {
        int(key): value for key, value in choices['exact_count_distribution'].items()}
    assert sum(row.get('reference_at_minimum') is True for row in rows) == choices['reference_at_minimum']
    assert sum(row.get('reference_at_minimum') is None for row in rows) == choices['reference_unknown']
    assert abs(np.mean([row['choices'] for row in exact])-choices['exact_subset_mean']) < 1e-12
    for row in exact:
        assert row['random_success'] == (1/row['choices'] if row['reference_at_minimum'] else 0)
    return paired, choices, rows


def build_comparisons(man, save):
    paired, choices, rows = comparison_evidence(man)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3), layout='constrained',
                             gridspec_kw={'width_ratios': [1.25, 1]})
    coverage, timing = axes
    policies = [('mappa', 'MAPPA sweep', '#087F8C'),
                ('slap_sweep', 'SLAP expanded sweep', '#CB702A'),
                ('slap_bidirectional', 'SLAP binary, uncut', '#586974')]
    for y, (key, label, color) in enumerate(policies):
        result = paired['coverage'][key]
        rate = 100*result['recovered']/paired['reactions']
        coverage.barh(y, rate, height=.4, color=color)
        coverage.text(2, y, f"{result['recovered']:,} / 1,851  ({rate:.2f}%)", va='center',
                      color='white', weight='bold', fontsize=9)
        coverage.text(1, y+.31, f"{result.get('not_recovered',0)} absent; "
                      f"{result.get('unknown',0)} unresolved", fontsize=8, color='#586974')
    coverage.set(yticks=range(3), yticklabels=[item[1] for item in policies], xlim=(0, 102),
                 ylim=(2.6, -.55), xlabel='Reference inclusion (%)')
    coverage.set_title('a  All 1,851 reactions · both directions', loc='left', weight='bold')
    stages = [('slap_standard', 'SLAP binary, uncut', '#586974'),
              ('mappa_search', 'MAPPA search', '#087F8C'),
              ('mappa_decode', 'MAPPA decoding', '#77559B'),
              ('slap_sweep', 'SLAP expanded sweep', '#CB702A')]
    for y, (key, label, color) in enumerate(stages):
        value = paired['paired_search'][key]['mean_cpu_seconds']
        timing.barh(y, value, color=color, height=.52)
        timing.text(value*1.12, y, f'{value:.3f}', va='center', fontsize=9)
    timing.set_xscale('log')
    timing.set(yticks=range(4), yticklabels=[item[1] for item in stages],
               xlim=(.03, 40), ylim=(3.6, -.55), xlabel='Mean CPU s / reaction (log scale)')
    timing.set_title('b  Same 1,834 completed reactions', loc='left', weight='bold')
    for ax in axes:
        ax.grid(axis='x', alpha=.15)
        ax.set_axisbelow(True)
    save(fig, 'fig9_paired_slap')

    exact = [row for row in rows if row['status'] == 'complete']
    bounds, labels = [1, 2, 5, 10, 20], ['1', '2', '3–5', '6–10', '11–20']
    groups = np.zeros((2, len(bounds)), dtype=int)
    for row in exact:
        column = next(i for i, bound in enumerate(bounds) if row['choices'] <= bound)
        groups[int(not row['reference_at_minimum']), column] += 1
    fig, ax = plt.subplots(figsize=(8.8, 4.2), layout='constrained')
    x = np.arange(len(labels))
    ax.bar(x, groups[0], color='#087F8C', label='Reference among choices')
    ax.bar(x, groups[1], bottom=groups[0], color='#CB702A', label='Reference absent')
    for i, total in enumerate(groups.sum(axis=0)):
        ax.text(i, total*1.15, f'{total:,}', ha='center', fontsize=9)
    ax.set_yscale('log')
    ax.set(xticks=x, xticklabels=labels, ylim=(.7, 5000),
           xlabel='Distinct minimum-event choices', ylabel='Reactions (log scale)')
    ax.set_title('MAPPA Golden · one seed · 1,829 certified counts', loc='left', weight='bold')
    ax.legend(frameon=False, loc='upper right', fontsize=9)
    ax.grid(axis='y', alpha=.15)
    ax.set_axisbelow(True)
    save(fig, 'fig10_minimum_choices')
