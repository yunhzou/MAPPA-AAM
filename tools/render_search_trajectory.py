"""Render recorded fragment-growth events in the shared 3D reaction style."""
import argparse
import hashlib
import json
from pathlib import Path
from mappa.presentation.trajectory import with_display_geometry
from mappa.presentation.viewers import growth_trace_html


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('trace',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--geometry-input',type=Path,
                   help='Original endpoint JSON with chemistry metadata for older traces')
    p.add_argument('--geometry-cache',type=Path,
                   help='Reuse or generate display conformers from indexed endpoint chemistry')
    a=p.parse_args()
    source=a.trace.read_bytes()
    document=json.loads(source)
    if a.geometry_input is not None or a.geometry_cache is not None:
        endpoints=None if a.geometry_input is None else json.loads(a.geometry_input.read_text())
        document=with_display_geometry(document,
            cache=a.geometry_cache or a.output.parent/'geometry',endpoints=endpoints)
        document['display_geometry']['source_trace_sha256']=hashlib.sha256(source).hexdigest()
        display_path=a.output.parent/'display-trace.json'
        if display_path.resolve()==a.trace.resolve():
            p.error('The raw trace must be separate from display-trace.json')
        display_path.write_text(json.dumps(document,separators=(',',':'))+'\n')
    a.output.write_text(growth_trace_html(document))
    print(a.output.resolve())
