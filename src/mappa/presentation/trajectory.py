"""Reusable, verified replay of fragment calls along saved AAM search paths.

Replay only saved fragment calls using the matcher's witness trace. No sweep
or seed-order search is run.
"""
import hashlib
import importlib
import gzip
import json
import math
from pathlib import Path

import numpy as np
from mappa.mappa_core.utils.io.checkpoints import read_aam_checkpoint
from mappa.mappa_core.utils.io.records import aam_from_record
from mappa.mappa_core.core.matching.graphs import build_graph
from mappa.mappa_core.post_processing.decoder.events import classify_bonds
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
growth = importlib.import_module('mappa.mappa_core.core.matching.growth.engine')


def replay(source, target, state, edge, po, ro, context):
    if context.growth_policy != 'growth_edge_splitmix64_v2':
        raise ValueError('Growth policy differs from archive; use its producing engine')
    events = []
    isos = growth.grow_island_python(source, target, edge.seed, dict(state.mapping),
        graph_floor=context.graph_floor, iso_tol=context.iso_tolerance,
        max_branches=context.branch_limit, events=events, random_seed=context.random_seed,
        islands_R=dict(state.islands), p_orbits=po, r_orbits=ro,
        prior_deferred_edges=state.deferred_edges)
    expected = dict(edge.match['symmetry']['witness'])
    matching = [i for i in isos if dict(i)==expected and sorted(i.fragment)==sorted(edge.match['fragment'])]
    if not matching or sorted(matching[0].deferred_edges)!=sorted(map(tuple,edge.match['deferred_edges'])):
        raise ValueError(f'Replay differs from archived transition {edge.id}; use its producing engine')
    return dict(events=events, isos=[dict(i) for i in isos],
                iso_candidates=[{**i.symmetry, 'symmetry_complete': True} for i in isos], verified=True)


def bond_events(raw, mapping, floor=.2, tolerance=.5, metal_tolerance=.3):
    """Same signed-event thresholds as decoding; floor only labels the drawing."""
    r=np.array(raw['reactant']['wbo']);p=np.array(raw['product']['wbo'])
    broken,formed,_,_=classify_bonds(mapping,r,p,dwbo_threshold=tolerance,
        elements_R=raw['reactant']['elements'],elements_P=raw['product']['elements'],
        metal_dwbo_threshold=metal_tolerance)
    inverse={p:r for r,p in mapping.items()}
    out=[]
    for a,b,x,y in broken:
        out.append(dict(kind='broken' if y is None or y<=floor else 'weakened',r=[a,b],
                        p=[mapping[r] for r in (a,b) if r in mapping],wbo=[x,y]))
    for u,v,x,y in formed:
        atoms=sorted(inverse[p] for p in (u,v) if p in inverse)
        targets=[mapping[r] for r in atoms] if len(atoms)==2 else [u,v]
        out.append(dict(kind='formed' if x is None or x<=floor else 'strengthened',r=atoms,
                        p=targets,wbo=[x,y]))
    return out


def frames_for(graph, terminal, traces):
    frames=[];stage=0
    path=next(graph.paths(terminal))
    atom_count=len(path.context.source_atoms)
    for stage,eid in enumerate(path.transitions,1):
        edge=graph.transitions[eid];state=graph.states[edge.source]
        if edge.match is None:
            after=graph.states[edge.target]
            frames.append(dict(kind='locked',stage=stage,title='Reuse recorded state',locked=dict(after.mapping),
                active=[],deferred=list(map(list,after.deferred_edges)),candidates=[],preferred=0))
            continue
        trace=traces[eid]
        locked=dict(state.mapping);deferred=list(map(list,state.deferred_edges))
        expected=dict(edge.match['symmetry']['witness'])
        patterns=[]
        for event in trace['events']:
            kind=event['type']
            if kind=='pop_skip':continue
            data=dict(stage=stage,seed=edge.seed,transition=eid,locked=locked,
                      deferred=[*deferred],kind=kind,expected=expected)
            if kind=='seed_start':
                patterns=event.get('cand_patterns',[])
                data.update(title=f'Start fragment {stage} at R{edge.seed}',active=event['fragment'],candidates=patterns)
            elif kind=='pop':
                trial_edge=[event['edge']['frag_atom'],event['edge']['ext_atom']]
                patterns=event['pre_state']['cands_pattern_sample']
                data.update(title=f'Test R{trial_edge[0]}–R{trial_edge[1]}',
                    active=event['pre_state']['fragment'],candidates=patterns,edge=trial_edge)
            elif kind in ('commit','consumed'):
                trial_edge=[event['edge']['frag_atom'],event['edge']['ext_atom']]
                patterns=event['cands_pattern_after'] if kind=='commit' else event['cands_pattern_sample']
                if kind=='consumed':deferred.append(trial_edge)
                data.update(title=(f'Grow to R{trial_edge[1]}' if kind=='commit' else
                    f'Defer R{trial_edge[0]}–R{trial_edge[1]}'),
                    active=event['fragment'],candidates=patterns,edge=trial_edge,
                    deferred=[*deferred])
            elif kind=='seed_end':
                patterns=trace['iso_candidates']
                data.update(title=f'Fragment {stage} saturated · {len(patterns)} retained placement(s)',
                    active=event['fragment'],candidates=patterns)
            else:continue
            data['candidate_kind']='compressed_match' if kind=='seed_end' else 'witness'
            # Growth frames show one witness; saturation shows every retained placement.
            data['preferred']=max(range(len(patterns)),key=lambda i:sum(expected.get(a)==b for a,b in patterns[i]['witness'].items())) if patterns else 0
            frames.append(data)
        after=graph.states[edge.target]
        frames.append(dict(kind='locked',stage=stage,seed=edge.seed,transition=eid,
            title=f'Commit fragment {stage} · {len(after.mapping)} / {atom_count} atoms assigned',
            locked=dict(after.mapping),active=[],deferred=list(map(list,after.deferred_edges)),
            candidates=[],preferred=0,expected=expected))
    final=dict(graph.states[terminal].mapping)
    frames.append(dict(kind='terminal',stage=stage,title='Saved terminal mapping · final bond events',
        locked=final,active=[],deferred=list(map(list,graph.states[terminal].deferred_edges)),
        candidates=[],preferred=0))
    return frames


def endpoint_record(endpoint):
    if endpoint.coordinates is None:
        raise ValueError('A 3D search trajectory requires endpoint coordinates in the archive')
    return dict(elements=list(endpoint.elements),coordinates=endpoint.coordinates.tolist(),wbo=endpoint.wbo.tolist())


def read_archive(path):
    if path.name.endswith(('.json','.json.gz')):
        opener=gzip.open if path.suffix=='.gz' else open
        with opener(path,'rt') as stream:return aam_from_record(json.load(stream))
    return read_aam_checkpoint(path)


def build_trajectory(selections, *, title=None, key_atoms=None, watch_targets=None,
                     event_tolerance=.5, metal_event_tolerance=.3):
    """Capture selected contexts from compatible archives of the same endpoints.

    Each selection specifies either an in-memory ``aam`` or an ``archive``,
    plus a context, with optional terminals, label,
    path_labels, default_terminal and focus={source_edge:[a,b], target:p}.
    One recorded history is followed per terminal; compressed permutations are
    represented, not expanded. A replay mismatch raises instead of publishing.
    """
    if not math.isfinite(event_tolerance) or event_tolerance <= 0:
        raise ValueError('event_tolerance must be finite and positive')
    if metal_event_tolerance is not None and (not math.isfinite(metal_event_tolerance) or metal_event_tolerance <= 0):
        raise ValueError('metal_event_tolerance must be finite and positive')
    runs,checks,raw=[],[],None
    for selection in selections:
        if ('aam' in selection) == ('archive' in selection):
            raise ValueError('Select exactly one of aam or archive')
        if 'aam' in selection:
            result=selection['aam']; origin='in_memory'
            provenance=dict(source='in_memory',problem=result.problem.name)
        else:
            archive=Path(selection['archive']);result=read_archive(archive);origin=str(archive)
            provenance=dict(archive=str(archive),sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
        g=result.graph
        context_id=int(selection.get('context',0));c=g.contexts[context_id]
        current=dict(name=result.problem.name,
                     reactant=endpoint_record(result.problem.reactant),product=endpoint_record(result.problem.product))
        if raw is None:raw=current
        elif any(raw[k]!=current[k] for k in ('reactant','product')):
            raise ValueError('Comparison archives must have identical endpoint atom order, coordinates and WBOs')
        source,target=[build_graph(raw[k]['elements'],np.array(raw[k]['wbo']),bond_cut=c.graph_floor)
                       for k in ('reactant','product')]
        source.remove_edges_from(c.cuts)
        po=_nauty_orbits(target,wbo_tol=c.iso_tolerance);ro=_nauty_orbits(source,wbo_tol=c.iso_tolerance)
        terminals=selection.get('terminals',[t for t in g.terminals if g.states[t].context==context_id])
        if not terminals:raise ValueError(f'Context {context_id} has no selected terminal paths')
        if any(g.states[t].context!=context_id or t not in g.terminals for t in terminals):
            raise ValueError('Selected terminal does not belong to the selected context')
        ids={eid for t in terminals for eid in next(g.paths(t)).transitions}
        traces={eid:replay(source,target,g.states[g.transitions[eid].source],g.transitions[eid],po,ro,c)
                for eid in sorted(ids) if g.transitions[eid].match is not None}
        paths=[]
        for t in terminals:
            frames=frames_for(g,t,traces);mapping=dict(g.states[t].mapping)
            complete=len(mapping)==len(raw['reactant']['elements'])
            events=bond_events(raw,mapping,c.graph_floor,event_tolerance,metal_event_tolerance) if complete else []
            for frame in frames:
                for candidate in frame['candidates']:
                    preview={**frame['locked'],**candidate['witness']}
                    if len(set(preview.values()))!=len(preview) or any(
                        raw['reactant']['elements'][a]!=raw['product']['elements'][b] for a,b in preview.items()):
                        raise ValueError('Invalid replayed candidate mapping')
            paths.append(dict(terminal=t,label=selection.get('path_labels',{}).get(str(t),f'Terminal {t}'),
                              frames=frames,mapping=mapping,events=events,complete=complete))
            checks.append(dict(archive=origin,context=context_id,terminal=t,frames=len(frames),
                               archived_fragment_matches='exact',complete=complete,events=len(events)))
        nodes=[s for s in g.states if s.context==context_id]
        runs.append(dict(label=selection.get('label',f'Context {context_id} · R{c.seed_order[0]} first'),
            context=context_id,seed_order=list(c.seed_order),cuts=list(c.cuts),paths=paths,
            default_terminal=selection.get('default_terminal',terminals[0]),focus=selection.get('focus'),
            config=dict(graph_floor=c.graph_floor,iso_tolerance=c.iso_tolerance,branch_limit=c.branch_limit,
                        growth_policy=getattr(c,"growth_policy","legacy_wbo_priority"),
                        random_seed=getattr(c,"random_seed",42),
                        event_tolerance=event_tolerance,metal_event_tolerance=metal_event_tolerance),
            graph=dict(states=[dict(id=s.id,assigned=len(s.mapping)) for s in nodes],
                transitions=[dict(id=eid,source=g.transitions[eid].source,target=g.transitions[eid].target,
                    seed=g.transitions[eid].seed,size=len(g.transitions[eid].match['fragment']) if g.transitions[eid].match else 0)
                    for eid in sorted(ids)]),
            provenance=provenance,traces=traces))
    if raw is None:raise ValueError('At least one archive selection is required')
    default_keys=sorted({a for run in runs for path in run['paths'] for event in path['events'] for a in event['r']})
    return dict(schema='rxn_core.search_trajectory/v1',name=title or raw['name'] or 'AAM search trajectory',
        input=raw,runs=runs,key_atoms=default_keys if key_atoms is None else key_atoms,
        watch_targets=[] if watch_targets is None else watch_targets,
        capture_source={str(p.relative_to(Path(__file__).parents[2])):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(__file__).parents[1].rglob('*.py'))},
        scope='Diagnostic replay of selected saved fragment calls; one recorded history per terminal. '
              'Growth frames show one feasible witness; saturation shows all completed compressed placements. '
              'Every fragment placement is checked against its archive. '
              'The endpoint coordinates stay fixed throughout playback.',checks=checks)


from copy import deepcopy


import json


from mappa.mappa_core.core.matching.graphs import build_graph


from mappa.mappa_core.core.matching.representation.state import candidate_from_record


from mappa.mappa_core.core.matching.symmetry.workspace import SymmetryWorkspace


def with_represented_targets(document):
    """Annotate per-atom marginal domains without enumerating full mappings.

    Local pools act first, then the conditioned candidate group, as in the
    family decoder. Marginal domains must never be shown as independently
    selectable assignments. Raw automorph-block unions are not target domains.
    """
    import numpy as np
    doc = deepcopy(document)
    ep = doc['input']['product']
    cache = {}
    for run in doc['runs']:
        cfg = run['config']
        target = build_graph(ep['elements'], np.asarray(ep['wbo']), bond_cut=cfg['graph_floor'])
        workspace = SymmetryWorkspace(target, cfg['iso_tolerance'])
        for path in run['paths']:
            for frame in path['frames']:
                for record in frame['candidates']:
                    if ((record.get('automorph_domains', 0) and 'automorph_blocks' not in record)
                            or (frame['kind'] == 'seed_end' and not record.get('symmetry_complete'))):
                        record['target_options'] = None
                        continue
                    key = json.dumps([cfg['graph_floor'], cfg['iso_tolerance'], frame['locked'], record], sort_keys=True)
                    if key not in cache:
                        cand = candidate_from_record(record)
                        locked = tuple(sorted((int(r), int(p)) for r, p in frame['locked'].items()))
                        generators = workspace.conditioned_generators(cand, locked)
                        block_r = {r for b in cand.blocks for r in b.r_atoms}
                        fixed = {p for _, p in locked} | {
                            p for r, p in cand.mapping.items() if r not in block_r or r in cand.exact_fixed}
                        options = {}
                        for r, p in cand.mapping.items():
                            domain = {p}
                            for block in cand.blocks:
                                free = set(block.p_atoms) - fixed
                                if domain & free:
                                    domain |= free
                            queue = list(domain)
                            while queue:
                                atom = queue.pop()
                                for gen in generators:
                                    image = gen[atom]
                                    if image not in domain:
                                        domain.add(image)
                                        queue.append(image)
                            options[str(r)] = sorted(domain)
                        cache[key] = options
                    record['target_options'] = cache[key]
    return doc

