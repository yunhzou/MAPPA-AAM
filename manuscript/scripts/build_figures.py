"""Build the final paper's vector figures from checked-in evidence; no AAM runs."""
import os
from pathlib import Path
MAN=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(MAN/'build/mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':240})
from build_molecule_figure import build as build_molecule_figure, build_design_examples, build_completion_example
build_molecule_figure(MAN)
build_design_examples(MAN)
build_completion_example(MAN)
from build_alternatives_figure import build as build_alternatives_figure
build_alternatives_figure(MAN)

from build_benchmark_figures import build as build_benchmark_figures
build_benchmark_figures(MAN)
print('Built current manuscript figures and tables from saved evidence.')
