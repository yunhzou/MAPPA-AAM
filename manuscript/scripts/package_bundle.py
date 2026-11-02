"""Bundle only material used by the final preprint, with a file-hash manifest."""
from pathlib import Path
import hashlib,json,zipfile
from publication_privacy import check_publication
MAN=Path(__file__).resolve().parents[1]
names=['manuscript.pdf','preprint.tex','preprint.bbl','references.bib','README.md','REPRODUCIBILITY.md','EDITORIAL_NOTES.md','requirements-build.txt','Makefile','natbib.sty','algorithmic.sty']
names += ['assets/'+n for n in ['matterlab.cls','paper-preamble.tex','plainnat.bst','Optimistic.ttf','logo_matterlab.pdf','logo_ac.pdf','logo_uoft.pdf','logo_vector.png','nvidia-logo-vert.png']]
names += ['includes/'+n for n in ['include-abstract.tex','include-acknowledgement.tex','include-body.tex','include-appendix.tex','paper.tex','supplement.tex','generated-cap-ablation-table.tex','generated-decoder-table.tex','generated-results.tex','generated-competitor-table.tex','generated-miss-table.tex','generated-direction-table.tex','generated-coordinate-minima-table.tex','generated-coordinate-time-table.tex','generated-end-to-end-table.tex','generated-end-to-end-text.tex','generated-end-to-end-scope.tex']]
names += ['figs/'+stem+'.'+ext for stem in ['fig1_algorithm','fig2_golden','fig3_coordinate','fig4_alternatives','fig5_pipeline','fig6_current_golden'] for ext in ['pdf','svg','png']]
names += ['evidence/'+n for n in ['current_golden.json','current_golden_cases.json','author-affiliations.md','paper_sources.json','seed_comparison.json','slap_sweep.json','final_dedup.json','PAPER_EVIDENCE.md','molecule_example.json','competitors.json','unswept.json','checkpoint_verdict_audit.json','golden_miss_analysis.json','direction_recovery.json','timing_comparison.json','coordinate_minima.json','golden_cap_ablation.json','multicandidate_example.json','output_multiplicity.json','golden_case9.json','golden_case9_verification.json','end_to_end_timing.json']]
names += ['scripts/'+n for n in ['build.sh','build_figures.py','validate_outputs.py','package_bundle.py','build_molecule_figure.py','build_alternatives_figure.py','publication_privacy.py']]
for n in names:check_publication(n,(MAN/n).read_bytes())
checks={n:hashlib.sha256((MAN/n).read_bytes()).hexdigest() for n in sorted(names)}
with zipfile.ZipFile(MAN/'manuscript_bundle.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for n in names:z.write(MAN/n,'manuscript/'+n)
 z.writestr('manuscript/SHA256SUMS.json',json.dumps(checks,indent=2)+'\n')
 for n in ['artifact-validation.json','visual-review.json']:
  if (MAN/'build'/n).is_file():
   record=json.loads((MAN/'build'/n).read_text())
   if record.get('manuscript_sha256')==checks['manuscript.pdf']:z.write(MAN/'build'/n,'manuscript/validation/'+n)
with zipfile.ZipFile(MAN/'manuscript_bundle.zip') as z:
 assert z.testzip() is None
 for name in z.namelist():check_publication(name,z.read(name))
print(f'Bundle: {len(names)} files, {(MAN/"manuscript_bundle.zip").stat().st_size/1e6:.2f} MB')
