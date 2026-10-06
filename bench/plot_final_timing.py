"""Plot recorded final MAPPA timings; no searches or decoding are run."""
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mappa-timing-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--report', type=Path,
    default=Path(__file__).resolve().parents[1] / 'reports/final_end_to_end_20260916')
args = parser.parse_args()
source = args.report / 'summary.json'
data = json.loads(source.read_text())
rows = sorted(data['rows'], key=lambda r: r['execution']['case'])
finished = [r for r in rows if r['execution']['complete']]
stopped = [r for r in rows if not r['execution']['complete']]
assert len(rows) == 140 and len(finished) == 138
assert {r['execution']['case'] for r in stopped} == {25, 31}
cpu = np.array([r['execution']['sampled_process_cpu_seconds'] for r in rows])
decode = np.array([r['result']['stages']['decode']['cpu_seconds'] for r in finished])
bins = np.histogram(decode, [0, 1, 10, 60, np.inf])[0]
assert sum(bins) == len(finished)
st = data['complete_cohort']
mean_total = st['end_to_end']['cpu_seconds']['mean']
means = [st[k]['cpu_seconds']['mean'] for k in ('search', 'catalogue', 'decode')]
means.append(mean_total - sum(means))
assert min(means) >= 0
ink, muted, blue, purple, red = '#193340', '#5A6C77', '#247EA2', '#8061AB', '#C3453C'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.labelcolor': ink, 'text.color': ink, 'xtick.color': muted,
    'ytick.color': muted, 'svg.fonttype': 'none'})
fig = plt.figure(figsize=(14.8, 10.6), facecolor='white')
grid = fig.add_gridspec(2, 2, left=.073, right=.968, top=.80, bottom=.125,
                       height_ratios=[1.05, 1], hspace=.55, wspace=.42)
fig.text(.073, .953, 'Where MAPPA spends its time', fontsize=25, weight='bold')
fig.text(.073, .916, 'Final version · 140 coordinate reactions · one-seed cut sweep · cap 2,000 · competition off',
         fontsize=11.5, color=muted)
fig.text(.073, .866, f'ALL ATTEMPTS   ≈{cpu.sum()/60:.1f} CPU min', fontsize=15, weight='bold', color=blue)
fig.text(.51, .866, '138 completed     |     2 stopped at the 5-minute watchdog', fontsize=12, color=red)

a = fig.add_subplot(grid[0, :])
for row in rows:
    e = row['execution']; x = e['case']; y = e['sampled_process_cpu_seconds']
    a.vlines(x, .07, y, color=red if not e['complete'] else '#D7E5EB', linewidth=1)
done_x = [r['execution']['case'] for r in finished]
done_y = [r['execution']['sampled_process_cpu_seconds'] for r in finished]
a.scatter(done_x, done_y, color=blue, s=23, zorder=3, label='Completed process')
a.scatter([r['execution']['case'] for r in stopped],
          [r['execution']['sampled_process_cpu_seconds'] for r in stopped],
          facecolors='white', edgecolors=red, marker='^', linewidths=1.8,
          s=100, zorder=5, label='Stopped: completion cost remains unknown')
for case, xtext, align in [(25, 21, 'right'), (31, 36, 'left')]:
    a.annotate(f'Case {case}', (case, cpu[case]), (xtext, 430),
               fontsize=10, color=red, ha=align,
               arrowprops={'arrowstyle': '-', 'color': red, 'linewidth': .8})
a.set(yscale='log', ylim=(.07, 650), xlim=(-2, 142), xlabel='Reaction index (all 140 attempts)',
      ylabel='Observed process CPU seconds\n(log scale)')
a.yaxis.set_major_locator(FixedLocator([.1, 1, 10, 100, 300]))
a.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
a.set_title('a   Full attempted cost, including interrupted runs', loc='left', fontsize=13, weight='bold', pad=12)
a.legend(loc='upper right', frameon=False, fontsize=9)
a.grid(axis='y', which='major', alpha=.17)
a.set_axisbelow(True)

a = fig.add_subplot(grid[1, 0])
positions = np.arange(4)
a.bar(positions, bins, color=[blue, '#5E9CB5', purple, red], width=.60)
for x, n in zip(positions, bins):
    a.text(x, n+2, str(n), ha='center', weight='bold', fontsize=12)
a.set(xticks=positions, xticklabels=['< 1 s', '1–10 s', '10–60 s', '≥ 60 s'],
      ylim=(0, 119), ylabel='Completed reactions', xlabel='Event-decoding CPU time')
a.set_title('b   Most decodes are fast; a few are costly', loc='left', fontsize=13, weight='bold', pad=36)
a.text(0, 1.04, f'138 completed only · median {np.median(decode):.3f} s · mean {np.mean(decode):.2f} s',
       transform=a.transAxes, fontsize=10, color=muted)
a.grid(axis='y', alpha=.15); a.set_axisbelow(True)

a = fig.add_subplot(grid[1, 1])
labels = ['Search', 'Family catalogue', 'Event decoding', 'Other / output']
a.barh(np.arange(4), means, color=[blue, '#6FAF9C', purple, '#ADB8BE'], height=.57)
for i, n in enumerate(means):
    a.text(n+.08, i, f'{n:.2f} s', va='center', fontsize=11,
           weight='bold' if i == 2 else 'normal', color=purple if i == 2 else ink)
a.set(yticks=np.arange(4), yticklabels=labels, xlim=(0, 6.7), ylim=(3.6, -.65),
      xlabel='Mean CPU seconds per completed reaction')
a.set_title('c   Decoding dominates completed-case cost', loc='left', fontsize=13, weight='bold', pad=36)
a.text(0, 1.04, f'138 completed only · end to end {mean_total:.2f} s mean · {st["end_to_end"]["cpu_seconds"]["median"]:.2f} s median',
       transform=a.transAxes, fontsize=10, color=muted)
a.grid(axis='x', alpha=.15); a.set_axisbelow(True)

fig.text(.073, .066, 'Scope: panel a samples whole-process CPU, including startup and stopped calls; its values are approximate.',
         fontsize=10, color=muted)
fig.text(.073, .043, 'Panels b–c measure the completed pipeline, excluding startup. No full-dataset completion mean is available: cases 25 and 31 are unfinished.',
         fontsize=10, color=muted)
for extension in ('png', 'svg'):
    fig.savefig(args.report / f'timing-overview.{extension}', dpi=180, facecolor='white')
svg = args.report / 'timing-overview.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
plt.close(fig)
metrics = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    attempts=len(rows), completed=len(finished), interrupted_cases=[25,31],
    sampled_process_cpu_seconds=float(cpu.sum()), sampled_process_cpu_minutes=float(cpu.sum()/60),
    sampled_process_cpu_mean_seconds=float(cpu.mean()),
    ten_costliest_share=float(np.sort(cpu)[-10:].sum()/cpu.sum()),
    completed_decode_bin_counts=dict(zip(['below_1s','1_to_10s','10_to_60s','at_least_60s'],map(int,bins))),
    completed_decode_cpu_share=means[2]/mean_total,
    distinction='Observed spending includes censored attempts; it is not the cost of completing all 140 decodes.')
(args.report / 'plot-metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')
print(json.dumps(metrics, indent=2))
