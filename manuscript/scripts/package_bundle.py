"""Export the single manuscript source as an Overleaf-ready project."""
from pathlib import Path
import zipfile

from publication_privacy import check_publication


MAN = Path(__file__).resolve().parents[1]
names = ['main.tex', 'references.bib', 'natbib.sty', 'algorithmic.sty', 'OVERLEAF.md']
names += ['assets/'+name for name in [
    'matterlab.cls', 'paper-preamble.tex', 'plainnat.bst', 'Optimistic.ttf',
    'logo_matterlab.pdf', 'logo_ac.pdf', 'logo_uoft.pdf', 'logo_vector.png',
    'nvidia-logo-vert.png',
]]
names += [str(path.relative_to(MAN)) for path in (MAN/'includes').glob('*.tex')]
names += [str(path.relative_to(MAN)) for path in (MAN/'figs').glob('*.pdf')]
names = sorted(set(names))
for name in names:
    check_publication(name, (MAN/name).read_bytes())

output = MAN/'mappa_overleaf.zip'
with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
    for name in names:
        archive.write(MAN/name, name)
with zipfile.ZipFile(output) as archive:
    assert archive.testzip() is None
    assert 'main.tex' in archive.namelist()
    for name in names:
        assert archive.read(name) == (MAN/name).read_bytes()
print(f'Overleaf: {len(names)} files, {output.stat().st_size/1e6:.2f} MB')
