"""Count saved fragment compositions, without symmetry queries or optimization.

Replay each fragment's recorded best witness in the saved action order. Count
its branch events directly from the original matrices and saved boundaries.
The result measures observed branches, not every attainable mapping or tie.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from mappa.mappa_core.utils.bonds import is_metal_element


def composition(pieces, mapping, reverse=False):
    pairs = [(tuple(sorted(atoms)), tuple(sorted(mapping[a] for a in atoms)))
             for atoms in pieces]
    return tuple(sorted((target, source) if reverse else (source, target)
                        for source, target in pairs))


def saved_branches(record, *, reverse=False, target_count=None):
    """Recount existing choices; never prepare a decoder or call a solver."""
    if record['schema'] != 'mappa.decoded_events/v14':
        raise ValueError('Require a saved binary-fragment decoded record')
    problem = record.get('event_problem', record['problem'])
    source = np.asarray(problem['reactant']['wbo'])
    target = np.asarray(problem['product']['wbo'])
    policy = record['event_policy']
    metals = np.array([is_metal_element(e) for e in problem['reactant']['elements']])
    thresholds = np.full(source.shape, policy['threshold'])
    if policy['metal_threshold'] is not None:
        thresholds[metals[:, None] | metals[None, :]] = policy['metal_threshold']
    catalogue = record['catalogue']
    prepared = {row['fragment']: row['optimization'] for row in record['fragment_responses']}
    operations = {}
    for fid, minimum in prepared.items():
        original = dict(catalogue['fragments'][fid]['witness'])
        witness = dict(minimum['witness'])
        operation = np.arange(len(target))
        for atom, image in original.items():
            operation[image] = witness[atom]
        if len(set(operation)) != len(operation):
            raise ValueError('Saved fragment witness is not an injective transport')
        operations[fid] = operation
    contexts = defaultdict(list)
    bounds = {tuple(key): value for key, value in record['window']['bounds']}
    for key in record['window']['contexts']:
        contexts[key[0]].append(tuple(key))
    counts, candidates = {}, set()
    evaluated = 0
    for saved in record['families']:
        fid, numbers = saved['family_id'], saved['tables']
        if len(numbers) != len(contexts[fid]):
            raise ValueError('Saved tables and retained contexts disagree')
        eligible = [(key, number) for key, number in zip(contexts[fid], numbers)
                    if target_count is None or bounds[key] <= target_count]
        if not eligible:
            continue
        family = catalogue['families'][fid]
        atoms, images = np.asarray(family['mapping'], dtype=int).T
        for kind, action in family['actions']:
            if kind == 'relation':
                images = operations[action][images]
            elif kind != 'group':
                raise ValueError('Unknown saved action')
            # An unchanged global placement is already a permitted saved choice.
        mapping = dict(zip(map(int, atoms), map(int, images)))
        if len(set(mapping.values())) != len(mapping):
            raise ValueError('Composed saved witness reuses a target atom')
        partitions = {}
        for context, number in eligible:
            table = record['tables'][number]
            pieces = tuple(tuple(part) for part in table['fragments'])
            if pieces not in partitions:
                pairs = [(a, b) for part in pieces for i, a in enumerate(part) for b in part[i+1:]]
                if pairs:
                    a, b = np.asarray(pairs, dtype=int).T
                    x, y = np.array([mapping[int(i)] for i in a]), np.array([mapping[int(i)] for i in b])
                    internal = np.count_nonzero(abs(source[a, b] - target[x, y]) > thresholds[a, b])
                else:
                    internal = 0
                partitions[pieces] = int(internal), composition(pieces, mapping, reverse)
            count, key = partitions[pieces]
            candidates.add(key)
            source_boundary = {tuple(map(int, pair)) for pair in table['reactant_boundary_pairs']}
            target_boundary = {tuple(map(int, pair)) for pair in table['product_boundary_pairs']}
            count += len(source_boundary) + len(target_boundary)
            owners = {a: i for i, part in enumerate(pieces) for a in part}
            for a, b in table['cut_bonds']:
                if a in owners and b in owners and owners[a] != owners[b]:
                    pair = tuple(sorted((mapping[a], mapping[b])))
                    count += int(abs(source[a, b] - target[pair]) > thresholds[a, b])
                    count -= int(tuple(sorted((a, b))) in source_boundary)
                    count -= int(pair in target_boundary)
            counts.setdefault(key, set()).add(count)
            evaluated += 1
    return counts, candidates, evaluated


def analyze(rows):
    from mappa.mappa_core.core.scheduler.bundles import Blob, read_blob
    records = [(row, read_blob(Blob(**row['decoded']), None).record) for row in rows]
    minimum = min(record['minimum_events'] for _, record in records)
    counts, candidates, contexts = {}, set(), 0
    for row, record in records:
        found, eligible, evaluated = saved_branches(record, reverse=row['reversed'], target_count=minimum)
        candidates.update(eligible)
        contexts += evaluated
        for key, values in found.items():
            counts.setdefault(key, set()).update(values)
    return dict(case=rows[0]['case'], detected_minimum=minimum,
        branches=sum(minimum in values for values in counts.values()),
        eligible_compositions=len(candidates), evaluated_contexts=contexts,
        observed_below_detected=sum(min(values) < minimum for values in counts.values()),
        decoded_hashes=[row['decoded']['sha256'] for row in rows])


def run(campaign, output, workers=1, cases=None):
    from decoded_reference import audit_inputs
    rows, hashes, _ = audit_inputs(campaign)
    grouped = defaultdict(list)
    for row in rows:
        if cases is None or row['case'] in cases:
            grouped[row['case']].append(row)
    if not grouped or any(sorted(r['reversed'] for r in pair) != [False, True] for pair in grouped.values()):
        raise ValueError('Require both directions for every reaction')
    started = time.perf_counter()
    report = dict(schema='mappa.saved_minimum_branches/v1', expected=len(grouped),
        original_hashes=hashes, runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        policy='Unique atom-set-pair compositions with a replayed saved local-witness composition at the saved detected minimum',
        additional_searches=0, results=[])
    output.parent.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for row in pool.map(analyze, grouped.values()):
            report['results'].append(row)
            if len(report['results']) % 100 == 0:
                print(json.dumps(dict(completed=len(report['results']), expected=len(grouped))), flush=True)
    report['results'].sort(key=lambda row: row['case'])
    report['elapsed_seconds'] = time.perf_counter()-started
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(completed=len(report['results']), zero_saved_branches=sum(r['branches']==0 for r in report['results']),
                         elapsed_seconds=report['elapsed_seconds'])), flush=True)


def publish(path, output):
    """Publish aggregate evidence and a figure, with no infrastructure fields."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report = json.loads(path.read_text())
    rows = report['results']
    if len(rows) != report['expected'] or len({row['case'] for row in rows}) != len(rows):
        raise ValueError('Require the complete requested reaction scope')
    values = np.array([row['branches'] for row in rows], dtype=int)
    summary = dict(reactions=len(rows), policy=report['policy'], additional_searches=report['additional_searches'],
        distribution=dict(sorted(Counter(map(int, values)).items())),
        mean=float(values.mean()), median=float(np.median(values)), p95=float(np.percentile(values,95)),
        maximum=int(values.max()), single_branch=int(np.count_nonzero(values==1)),
        no_saved_witness_at_minimum=int(np.count_nonzero(values==0)),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), runner_sha256=report['runner_sha256'])
    output.mkdir(parents=True, exist_ok=True)
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with gzip.open(output/'reactions.json.gz','wt') as stream:
        json.dump(rows,stream)
    edges = [0,1,2,5,10,20,50,100,500,float('inf')]
    labels = ['0','1','2','3–5','6–10','11–20','21–50','51–100','101–500','>500']
    bins = Counter(next(i for i,edge in enumerate(edges) if value<=edge) for value in values)
    used = [i for i in range(max(bins)+1) if i or bins[0]]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig, ax = plt.subplots(figsize=(8.6,4.5), layout='constrained')
    heights = [bins[i] for i in used]
    ax.bar(range(len(used)),heights,color='#087F8C',width=.7)
    for x,total in enumerate(heights):
        if total:
            ax.annotate(f'{total:,}',(x,total),xytext=(0,5),textcoords='offset points',ha='center')
    ax.set(xticks=range(len(used)),xticklabels=[labels[i] for i in used],yscale='log',
        ylim=(.7,max(heights)*2.2),ylabel='Golden reactions (log scale)',
        xlabel='Distinct branches with a recorded witness at the detected minimum')
    ax.set_title(f'Golden · {len(rows):,} reactions · one seed',loc='left',weight='bold',pad=14)
    ax.text(1,1.02,'Internal symmetry arrangements count once',transform=ax.transAxes,
            ha='right',va='bottom',color='#586974',fontsize=9)
    ax.grid(axis='y',alpha=.15); ax.set_axisbelow(True)
    for suffix in ('png','pdf','svg'):
        fig.savefig(output/f'minimum-branches.{suffix}',dpi=220,facecolor='white')
    plt.close(fig)
    print(json.dumps(summary))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--cases', type=int, nargs='+')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.plot:
        publish(args.campaign,args.output)
    else:
        run(args.campaign, args.output, args.workers, args.cases)


if __name__ == '__main__':
    main()
