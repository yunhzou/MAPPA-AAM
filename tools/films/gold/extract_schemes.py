"""Extract unaltered pathway-scheme regions from the supplied reference PDF.

Only the selected scheme regions (including original titles and footnotes)
are rasterized. The paper itself is not copied into the repository.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

REGIONS = [
    dict(file='scheme-2-route-a', scheme=2, route='a', page=4, journal_page=2985,
         crop=[164, 174, 2096, 962]),
    dict(file='scheme-4-route-b', scheme=4, route='b', page=5, journal_page=2986,
         crop=[164, 1030, 2096, 940]),
    dict(file='scheme-5-route-c', scheme=5, route='c', page=6, journal_page=2987,
         crop=[164, 1300, 2096, 962]),
]


def extract(pdf, output):
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for region in REGIONS:
        x, y, width, height = region['crop']
        target = output / region['file']
        subprocess.run(['pdftoppm', '-f', str(region['page']), '-l', str(region['page']),
                        '-scale-to-x', '2408', '-scale-to-y', '3200',
                        '-x', str(x), '-y', str(y), '-W', str(width), '-H', str(height),
                        '-singlefile', '-png', str(pdf), str(target)],
                       check=True, capture_output=True, timeout=45)
        png = target.with_suffix('.png')
        records.append({**region, 'file': png.name,
                        'sha256': hashlib.sha256(png.read_bytes()).hexdigest()})
    manifest = dict(source_title='Mechanism of the Gold-Catalyzed Rearrangement of '
                    '(3-Acyloxyprop-1-ynyl)oxiranes: A Dual Role of the Catalyst',
                    doi='10.1021/jo802516k', source_pdf=pdf.name,
                    source_pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
                    method='Direct PDF-region rasterization; original drawings, titles, and footnotes retained.',
                    page_render_pixels=[2408, 3200], schemes=records)
    (output / 'provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdf', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    extract(args.pdf, args.output)
