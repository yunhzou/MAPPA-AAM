"""Shared molecular and batch presentation; search code does not depend on it.

The molecular layout is the original white R/P/TS viewer from stable/legacy-aam.
Adapters supply recorded data, rather than introducing report-specific skins.
"""
from mappa.mappa_core.domain import AAMResult
from mappa.mappa_core.utils.io.records import aam_json
from .interpolation import proper_align_coordinates

import html
import json
from pathlib import Path
import re
from functools import wraps

import numpy as np

ASSETS = Path(__file__).with_name('static')


def stylesheet(family='reaction'):
    if family not in ('reaction', 'batch'):
        raise ValueError('Supported viewer styles are reaction and batch')
    return (ASSETS / f'{family}_viewer.css').read_text()


def style_document(document, family='reaction', layout='standard'):
    """Apply one registered style to an existing structural layout."""
    document = re.sub(r'<style[^>]*>.*?</style>', '', document, flags=re.S)
    style = f'<style data-viewer-style="{family}">{stylesheet(family)}</style>'
    if '</head>' in document:
        document = document.replace('</head>', style + '</head>', 1)
    else:
        document = re.sub(r'<!doctype[^>]*>', '', document, flags=re.I)
        document = '<!doctype html>' + style + document
    if re.search(r'<body(?:\s[^>]*)?>', document):
        document = re.sub(r' data-viewer-layout="[^"]*"', '', document, count=1)
        document = re.sub(r'<body(?=[\s>])', f'<body data-viewer-layout="{layout}"', document, count=1)
    else:
        document = re.sub(r'(?=<(?:header|main|h1|h2|div|section|article)\b)',
                          '<body data-viewer-layout="' + layout + '">', document, count=1)
    return document


def _json(value):
    return json.dumps(value, separators=(',', ':'), default=lambda x: x.tolist()).replace('</', '<\\/')


def viewer_style(family='reaction', layout='standard'):
    def decorate(function):
        @wraps(function)
        def render(*args, **kwargs):
            return style_document(function(*args, **kwargs), family, layout)
        return render
    return decorate


def align_product(reactant, product):
    """Rigid display fitting only; retain the recorded atom correspondence."""
    return proper_align_coordinates(product, reactant).tolist()


def comparison_document(case):
    """Adapt complete saved R/P or endpoint/TS witnesses to the original viewer."""
    r, p = case['endpoints']
    mechanisms = []
    for index, record in enumerate(case['records'], 1):
        source = record['mapping']
        mapping = dict(source) if source and isinstance(source[0], (tuple, list)) else dict(enumerate(source))
        assert sorted(mapping) == sorted(mapping.values()) == list(range(len(r['elements'])))
        assert all(r['elements'][a] == p['elements'][b] for a, b in mapping.items())
        ordered = [p['coordinates'][mapping[a]] for a in range(len(mapping))]
        # Decoder signs describe WBO decreases/increases, including bonds that
        # remain present. Presentation uses the same >0.2 bond-presence cutoff
        # as the endpoint renderer to distinguish loss/gain from order changes.
        events = []
        for event in record.get('events', []):
            event = dict(event)
            if event['kind'] in ('broken', 'formed') and 'wbo' in event:
                before, after = event['wbo']
                if before is not None and after is not None and before > 0.2 and after > 0.2:
                    event['decoder_kind'] = event['kind']
                    event['kind'] = 'weakened' if before > after else 'strengthened'
            events.append(event)
        by_kind = {kind: [e['r'] for e in events if e['kind'] == kind]
                   for kind in ('broken', 'formed', 'strengthened', 'weakened', 'order_changed')}
        mechanisms.append(dict(
            id=index, label=record['name'], mapping_RP=mapping,
            product_xyz_in_R=ordered, product_xyz_in_R_aligned=align_product(r['coordinates'], ordered),
            broken_bonds_R=by_kind['broken'], formed_bonds_R=by_kind['formed'],
            formed_bonds_P=[[mapping[a], mapping[b]] for a, b in by_kind['formed']],
            order_events={k: by_kind[k] for k in ('strengthened', 'weakened', 'order_changed')},
            event_records=events, event_counts={kind: len(pairs) for kind, pairs in by_kind.items()},
            pattern_id=record.get('pattern'), provenance=record.get('provenance', {}),
            core_atoms=sorted({a for e in events for a in e['r']}), igs=[], gt=None,
        ))
    return dict(step=str(case['name']), index=case['index'], n_atoms=len(r['elements']),
                default_mech_id=1, include_gt=False, note=case.get('note', ''),
                reactant=dict(elements=r['elements'], coords=r['coordinates'], wbo=r['wbo']),
                product=dict(elements=p['elements'], coords=p['coordinates'], wbo=p['wbo']),
                mechanisms=mechanisms)


def reaction_html(document):
    template = (ASSETS / 'reaction_viewer.html').read_text()
    return (template.replace('__TITLE__', html.escape(document['step']))
            .replace('__VIEWER_STYLE__', stylesheet())
            .replace('__LIBRARY__', (ASSETS / '3Dmol-min.js').read_text())
            .replace('__ZIP_LIBRARY__', (ASSETS / 'jszip.min.js').read_text())
            .replace('__DATA__', _json(document)))


def aam_growth_html(aam, *, context=None, terminals=None, title=None,
                    event_threshold=0.5, metal_event_threshold=0.3):
    """Offline interactive growth replay directly from an AAMResult.

    ``context=None`` includes every saved sweep/seed context having terminals;
    an integer selects one. No terminal/history or symmetry enumeration cap is
    added. One recorded history per terminal is replayed against its producer's
    matcher and checked exactly; this is not exhaustive decoded-event output.
    Coordinates are required. Capture temporarily observes matcher functions:
    call serially, not concurrently with another search/replay in this process.
    HTML embeds its renderer and data. In Jupyter, display it in an iframe so
    several viewers can coexist without sharing JavaScript state.
    """
    from mappa.presentation.trajectory import build_trajectory
    if terminals is not None and context is None:
        raise ValueError('Select one context when specifying terminals')
    if context is None:
        contexts=sorted({aam.graph.states[t].context for t in aam.graph.terminals})
    else:
        if not isinstance(context,int) or not 0 <= context < len(aam.graph.contexts):
            raise ValueError('context must index a saved search context')
        contexts=[context]
    selections=[]
    for i in contexts:
        c=aam.graph.contexts[i]
        label=f'Sweep {i}: '+(', '.join(f'R{a}–R{b}' for a,b in c.cuts) if c.cuts else 'uncut')
        selection=dict(aam=aam,context=i,label=label)
        if terminals is not None:selection['terminals']=list(terminals)
        selections.append(selection)
    document=build_trajectory(selections,title=title,event_tolerance=event_threshold,
                              metal_event_tolerance=metal_event_threshold)
    return growth_trace_html(document)


def growth_trace_html(document):
    """Present recorded algorithm events using the shared 3D reaction style."""
    from mappa.presentation.trajectory import with_represented_targets
    document = with_represented_targets(document)
    template = (ASSETS / 'growth_trace.html').read_text()
    page = style_document(template, layout='growth_trace')
    return (page.replace('__LIBRARY__', (ASSETS / '3Dmol-min.js').read_text())
            .replace('__DATA__', _json(document)))


def collection_html(documents, title='R/P and TS comparison'):
    """One copy of the shared renderer and libraries, with a saved-case selector."""
    if not documents:
        raise ValueError('No cases to display')
    template = ((ASSETS / 'reaction_viewer.html').read_text()
                .replace('__VIEWER_STYLE__', stylesheet())
                .replace('__LIBRARY__', (ASSETS / '3Dmol-min.js').read_text())
                .replace('__ZIP_LIBRARY__', (ASSETS / 'jszip.min.js').read_text()))
    page = '''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__PAGE_TITLE__</title></head><body><header class="topbar"><b>__PAGE_TITLE__</b>
<label>Case <select id="case"></select></label></header><iframe id="reaction" title="Original R/P/TS viewer"></iframe>
<script>const CASES=__CASES__,TEMPLATE=__TEMPLATE__,select=document.getElementById('case'),frame=document.getElementById('reaction');
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
CASES.forEach((c,i)=>{const o=document.createElement('option');o.value=i;o.textContent=(c.index===undefined?'':c.index+' · ')+c.step;select.appendChild(o)});
function render(){const c=CASES[+select.value];frame.srcdoc=TEMPLATE.replaceAll('__TITLE__',escape(c.step)).replace('__DATA__',JSON.stringify(c).replaceAll('</','<\\/'));history.replaceState(null,'','#case='+encodeURIComponent(c.index===undefined?select.value:c.index));}
const query=new URLSearchParams(location.hash.slice(1)).get('case'),choice=CASES.findIndex(c=>String(c.index)===query);if(choice>=0)select.value=choice;
select.onchange=render;render();</script></body></html>'''
    page = page.replace('__PAGE_TITLE__', html.escape(title)).replace('__CASES__', _json(documents)).replace('__TEMPLATE__', _json(template))
    # Avoid treating the embedded template's encoded styles as outer styles.
    return page.replace('</head>', '<style data-viewer-style="reaction">' + stylesheet() + '</style></head>', 1).replace('<body>', '<body data-viewer-layout="collection">', 1)


def write_aam_bundle(result: AAMResult, output_directory):
    """Save raw AAM and an offline graph/path viewer; never rerun matching."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    serialized = aam_json(result)
    (output / 'aam.json').write_text(serialized)
    assets = Path(__file__).parent / 'static'
    page = (assets / 'aam_search.html').read_text()
    page = style_document(page, layout='search_graph')
    page = page.replace('__LIBRARY__', (assets / '3Dmol-min.js').read_text())
    page = page.replace('__DATA__', serialized.replace('<', '\\u003c'))
    (output / 'search.html').write_text(page)
    return output

