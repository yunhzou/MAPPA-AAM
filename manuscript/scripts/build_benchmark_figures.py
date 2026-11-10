"""Current benchmark figures and tables from publication-safe saved evidence."""
import gzip
import hashlib
import json
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter


def evidence(man):
    summary = json.loads((man/'evidence/current_benchmarks.json').read_text())
    offsets = json.loads(gzip.decompress((man/'evidence/reference_offsets.json.gz').read_bytes()))
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
        assert row['difference'] == row['reference_count']-row['overall_minimum']
        assert row['source'] in {'verified', 'anchored'}
        assert row['certificate_sha256'] or row['source'] == 'anchored'
        assert not {'branch_minimum', 'offset_lower', 'offset_upper'} & row.keys()
    assert offsets['additional_searches'] == 0
    return summary, offsets


def offset_counts(rows):
    """One bin per verified reference completion; no branch-minimum estimation."""
    return Counter(row['difference'] for row in rows)


def anchored_evidence(man):
    """Validate the saved completions of references absent from the original output."""
    data = json.loads((man/'evidence/anchored_offsets.json').read_text())
    rows = data['results']
    offsets = json.loads(gzip.decompress((man/'evidence/branch_offsets.json.gz').read_bytes()))
    absent = {(r['dataset'], r['case']) for r in offsets['results']
              if r['status'] == 'not_recovered'}
    assert len(rows) == len(absent) == 196
    assert {(r['dataset'], r['case']) for r in rows} == absent
    assert Counter(r['dataset'] for r in rows) == {'Golden': 18, 'EnzymeMap': 178}
    for row in rows:
        pair = row['directions']
        assert len(pair) == 2 and {r['reversed'] for r in pair} == {False, True}
        assert all(r['all_anchors_preserved'] and r['maximum_atom_usage'] for r in pair)
        assert all(r['model_event_count'] == r['physical_event_count'] for r in pair)
        assert row['anchored_count'] == min(r['model_event_count'] for r in pair)
        assert row['difference'] == row['anchored_count']-row['overall_minimum']
    return data


def update_reference_offsets(man, enzyme_counts):
    """Extract plotted counts from existing certificates, without querying MAPPA."""
    paths = dict(
        baseline=man/'evidence/branch_offsets.json.gz',
        golden=man/'evidence/decoded_verification_closure_directions.json.gz',
        anchored=man/'evidence/anchored_offsets.json',
        enzyme=enzyme_counts,
    )
    original = json.loads(gzip.decompress(paths['baseline'].read_bytes()))
    golden = json.loads(gzip.decompress(paths['golden'].read_bytes()))['studies']['seeds1']
    enzyme = json.loads(gzip.decompress(enzyme_counts.read_bytes()))
    anchored = anchored_evidence(man)
    summary = json.loads((man.parent/'reports/enzymemap_native_20261106/decoded-reference-validation.json').read_text())
    assert enzyme['records_sha256'] == summary['records_sha256']
    assert enzyme['positive_directions'] == summary['positive_certificates_checked']

    selected = {}
    for dataset, records in [('Golden', golden['directions']), ('EnzymeMap', enzyme['results'])]:
        for row in records:
            if dataset == 'Golden' and row['status'] != 'recovered':
                continue
            key = dataset, row['case']
            if key not in selected or (row['event_count'], row['reversed']) < (
                    selected[key]['event_count'], selected[key]['reversed']):
                selected[key] = row
    assert Counter(dataset for dataset, _ in selected) == {'Golden': 1833, 'EnzymeMap': 47796}
    completions = {(r['dataset'], r['case']): r for r in anchored['results']}
    results = []
    for old in original['results']:
        key = old['dataset'], old['case']
        if old['status'] == 'not_recovered':
            assert key not in selected
            completion = completions[key]
            row = min(completion['directions'], key=lambda r: (r['model_event_count'], r['reversed']))
            baseline = completion['overall_minimum']
            count = row['model_event_count']
            source, certificate, decoded, reference = 'anchored', None, None, None
        else:
            row = selected[key]
            baseline, count = old['overall_minimum'], row['event_count']
            source = 'verified'
            certificate, decoded, reference = (row[field] for field in
                ('certificate_sha256', 'decoded_sha256', 'reference_sha256'))
            assert isinstance(count, int) and len(certificate) == 64
        results.append(dict(dataset=key[0], case=key[1], source=source,
            overall_minimum=baseline, reference_count=count, difference=count-baseline,
            reversed=row['reversed'], certificate_sha256=certificate,
            decoded_sha256=decoded, reference_sha256=reference))

    data = dict(schema='mappa.reference_event_offsets/v1',
        metric='Verified GT-consistent mapping event count minus saved detected minimum',
        counting='Same saved all-atom independent-fragment event model on both sides',
        selection='Lowest count among saved verified forward/reverse witnesses; forward on ties; no new optimization',
        completion='Use saved anchored completions for the 196 absent references; independent coverage unchanged',
        additional_searches=0,
        source_sha256={name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()},
        enzyme_records_sha256=enzyme['records_sha256'], results=results)
    output = man/'evidence/reference_offsets.json.gz'
    output.write_bytes(gzip.compress(json.dumps(data, separators=(',', ':')).encode(), mtime=0))
    return data


def build_reference_offsets(man, offsets, save):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.9), layout='constrained')
    for column, (name, color) in enumerate(zip(('Golden', 'EnzymeMap'), ('#087F8C', '#3366A5'))):
        cohort = [r for r in offsets['results'] if r['dataset'] == name]
        counts = offset_counts(cohort)
        assert sum(counts.values()) == len(cohort)
        ax = axes[column]
        xs = sorted(counts)
        ax.bar(xs, [counts[x] for x in xs], width=.68, color=color)
        for x in xs:
            ax.text(x, counts[x]*1.2, f'{counts[x]:,}', ha='center', va='bottom', fontsize=8)
        first, last = min(0, min(xs)), max(xs)
        ax.set(xlim=(first-.7, last+.7), ylim=(.7, max(counts.values())*5),
               xlabel='Event offset n from detected minimum',
               xticks=[x for x in range(first, last+1) if x >= 0 or x in counts])
        ax.set_title(f"{'ab'[column]}  {name} · {len(cohort):,} reactions",
                     loc='left', weight='bold', pad=12)
    for ax in axes.flat:
        ax.set_yscale('log')
        ax.set_ylabel('Reactions (log scale)')
        ax.axvline(0, color='#223340', linewidth=.7, alpha=.5, zorder=0)
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    save(fig, 'fig7_reference_offsets')


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
        found = study['catalogue']['recovered']
        percent = 100*found/study['reactions']
        coverage.barh(y, percent, height=.36, color=colors[0])
        coverage.text(2, y, f"{found:,} / {study['reactions']:,}  ({percent:.3f}%)",
                      color='white', weight='bold', va='center', fontsize=9)
        coverage.text(1, y+.29, f"{study['catalogue']['not_recovered']:,} not recovered", fontsize=8, color='#586974')
    coverage.set(yticks=range(3), yticklabels=labels, xlim=(0, 102), ylim=(2.65, -.5),
                 xlabel='Confirmed reactions (%)', xticks=[0, 25, 50, 75, 100])
    coverage.set_title('a  Full-output reference coverage', loc='left', weight='bold', pad=12)
    left = np.zeros(3)
    stages = [(('matching', 'deduplication'), 'Matching', colors[0]),
              (('decoding',), 'Decoding', colors[2])]
    for keys, label, color in stages:
        values = np.array([sum(s['cpu_seconds'][key] for key in keys)/s['reactions']
                           for s in studies])
        cost.barh(range(3), values, left=left, height=.48, color=color, label=label)
        left += values
    for y, value in enumerate(left):
        cost.text(value+.25, y, f'{value:.2f}', va='center', fontsize=9)
    cost.set(yticks=range(3), yticklabels=labels, xlim=(0, max(left)*1.18), ylim=(2.65, -.8),
             xlabel='Mean compute CPU seconds / reaction')
    cost.set_title('b  Matching and decoding', loc='left', weight='bold', pad=12)
    cost.legend(loc='upper left', ncol=2, fontsize=8, frameon=False)
    for ax in (coverage, cost):
        ax.grid(axis='x', alpha=.15)
        ax.set_axisbelow(True)
    save(fig, 'fig6_current_benchmarks')

    build_reference_offsets(man, offsets, save)

    assert all(s['decoded']['unknown'] == 0 for s in studies)
    coverage_rows, timing_rows = [], []
    for study in studies:
        name = f"{study['dataset']}, {study['seeds']} seed"+('s' if study['seeds'] > 1 else '')
        def rate(count):
            return f"{count:,} ({100*count/study['reactions']:.3f}\\%)"
        coverage_rows.append(' & '.join([name, f"{study['reactions']:,}",
            rate(study['catalogue']['recovered']), rate(study['decoded']['recovered']),
            str(study['decoded']['not_recovered'])])+r'\\')
        matching, reduction, decoding = [study['cpu_seconds'][s]/study['reactions']
                                         for s in ('matching', 'deduplication', 'decoding')]
        timing_rows.append(' & '.join([name, *[f'{v:.2f}' for v in
            (matching+reduction, decoding, matching+reduction+decoding)]])+r'\\')
    def table(name, columns, header, rows):
        text = r'\begin{tabular}{'+columns+r'}\toprule'+'\n'+header+r'\\\midrule'+'\n'
        text += '\n'.join(rows)+'\n'+r'\bottomrule\end{tabular}'+'\n'
        (man/'includes'/name).write_text(text)
    table('generated-current-coverage-table.tex', 'lrrrr',
          r'Study & Reactions & Catalogue & Decoded & \shortstack{Decoded\\absent}', coverage_rows)
    table('generated-current-timing-table.tex', 'lrrr',
          'Study & Matching & Decoding & Total', timing_rows)
    paired, _, _ = comparison_evidence(man)
    competitors = json.loads((man/'evidence/competitors.json').read_text())
    method_rows = []
    def method_row(name, setting, output, found):
        method_rows.append(' & '.join([name, setting, output,
            f'{found:,} ({100*found/1851:.2f}\\%)'])+r'\\')
    for study in studies:
        if study['dataset'] == 'Golden':
            setting = f"{study['seeds']} seed"+('s' if study['seeds'] > 1 else '')+', sweep'
            method_row('MAPPA', setting, 'Families', study['catalogue']['recovered'])
    for key, setting in [('slap_sweep', 'Expanded sweep'),
                         ('slap_bidirectional', 'Binary, both directions')]:
        result = paired['coverage'][key]
        method_row('SLAP', setting, 'Candidates', result['recovered'])
    names = {'slap_binary': 'SLAP', 'localmapper': 'LocalMapper', 'chython': 'Chython',
             'rxnmapper': 'RXNMapper', 'rdt': 'RDT', 'indigo': 'Indigo'}
    for key, name in names.items():
        result = next(row for row in competitors['methods'] if row['method'] == key)
        method_row(name, 'Binary, default' if key == 'slap_binary' else 'Default',
                   'Candidates' if key == 'slap_binary' else 'Single mapping',
                   result['any_correct'])
    table('generated-method-coverage-table.tex', '@{}lllr@{}',
          'Method & Search setting & Output & Coverage', method_rows)
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
