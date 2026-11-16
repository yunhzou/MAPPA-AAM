"""Render saved comparison JSON through the shared original R/P/TS presentation."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import webbrowser

from mappa.presentation.viewers import collection_html, comparison_document


def render(source, output, *, golden_cases=None, geometry_cache=None):
    output = Path(output)
    if golden_cases is not None:
        from golden import reference_view_cases
        cases = reference_view_cases(source, golden_cases,
            geometry_cache=geometry_cache or output.parent/'geometry')
    else:
        cases = json.loads(Path(source).read_text())
    documents = []
    for case in cases:
        documents.append(comparison_document(case))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(collection_html(documents))
    if golden_cases is not None:
        output.with_suffix('.json').write_text(json.dumps(cases, indent=2)+'\n')
    return len(documents)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--golden-cases', nargs='+', type=int,
                        help='Render selected GT cases from a saved Golden input index')
    parser.add_argument('--geometry-cache', type=Path)
    parser.add_argument('--open', action='store_true', help='Open the finished HTML in the browser')
    args = parser.parse_args()
    started = time.perf_counter()
    count = render(args.source, args.output, golden_cases=args.golden_cases,
                   geometry_cache=args.geometry_cache)
    print(f'Rendered {count} saved cases in {time.perf_counter()-started:.3f} s; no matching or decoding.')
    print(args.output.resolve())
    if args.open:
        if sys.platform == 'darwin':
            subprocess.run(['open', args.output.resolve().as_uri()], check=True)
        elif not webbrowser.open(args.output.resolve().as_uri()):
            raise RuntimeError('The browser could not open the saved viewer')
