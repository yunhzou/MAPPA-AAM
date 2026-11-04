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
growth = importlib.import_module('mappa.mappa_core.core.matching.growth.engine')


def _same_placement_group(left, right, degree):
    if left == right:
        return True
    from sympy.combinatorics import Permutation, PermutationGroup
    identity = Permutation(list(range(degree)))
    groups = [PermutationGroup([Permutation(list(g)) for g in generators] or [identity])
              for generators in (left, right)]
    return all(group.contains(generator) for group, other in
               ((groups[0], groups[1]), (groups[1], groups[0]))
               for generator in other.generators)


def replay(source, target, state, after, edge, context):
    if context.growth_policy != 'growth_edge_splitmix64_v2':
        raise ValueError('Growth policy differs from archive; use its producing engine')
    events = []
    result = growth.grow_island_python(source, target, edge.seed, dict(state.mapping),
        graph_floor=context.graph_floor, iso_tol=context.iso_tolerance,
        max_branches=context.branch_limit, events=events, random_seed=context.random_seed,
        islands_R=dict(state.islands), self_tolerance=context.self_tolerance,
        prior_deferred_edges=state.deferred_edges)
    isos = result.matches
    matching = [fragment for fragment in isos if fragment.key == edge.match.key]
    if (not matching or sorted(result.deferred_edges) != sorted(after.deferred_edges)
            or not _same_placement_group(result.placement_generators,
                                         edge.placement_generators, len(target))):
        raise ValueError(f'Replay differs from archived transition {edge.id}; use its producing engine')
    return dict(events=events, isos=[dict(fragment.witness) for fragment in isos],
                placement_generators=result.placement_generators,
                iso_candidates=[dict(witness=dict(fragment.witness), blocks=[],
                    relation=fragment.to_record()) for fragment in isos], verified=True)


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
        expected=dict(edge.match.witness)
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
                data.update(title=f'Fragment {stage} saturated · {len(patterns)} canonical placement(s)',
                    active=event['fragment'],candidates=patterns,
                    placement_generators=trace['placement_generators'])
            else:continue
            data['candidate_kind']='compressed_match' if kind=='seed_end' else 'witness'
            # Saturation retains local relations and their shared placement action.
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
    return dict(elements=list(endpoint.elements),coordinates=endpoint.coordinates.tolist(),
                wbo=endpoint.wbo.tolist(),metadata=dict(endpoint.metadata))


def with_display_geometry(document, *, cache, endpoints=None):
    """Replace coordinates only in a rendering copy; retain the scientific trace."""
    from mappa.presentation.geometry import display_endpoint

    endpoints = document['input'] if endpoints is None else endpoints
    displayed = dict(document, input=dict(document['input']))
    metadata = {}
    for side in ('reactant', 'product'):
        original, endpoint = document['input'][side], endpoints[side]
        if (tuple(endpoint['elements']) != tuple(original['elements'])
                or not np.array_equal(endpoint['wbo'], original['wbo'])):
            raise ValueError('Display geometry input must preserve endpoint atom order and WBOs')
        rendered = display_endpoint(endpoint, cache=cache)
        displayed['input'][side] = dict(original, coordinates=rendered['coordinates'])
        metadata[side] = rendered['geometry']
    displayed['display_geometry'] = dict(endpoints=metadata,
        source_input_sha256=hashlib.sha256(json.dumps(document['input'], sort_keys=True,
            separators=(',', ':')).encode()).hexdigest())
    return displayed


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
        terminals=selection.get('terminals',[t for t in g.terminals if g.states[t].context==context_id])
        if not terminals:raise ValueError(f'Context {context_id} has no selected terminal paths')
        if any(g.states[t].context!=context_id or t not in g.terminals for t in terminals):
            raise ValueError('Selected terminal does not belong to the selected context')
        ids={eid for t in terminals for eid in next(g.paths(t)).transitions}
        traces={eid:replay(source,target,g.states[g.transitions[eid].source],
                          g.states[g.transitions[eid].target],g.transitions[eid],c)
                for eid in sorted(ids) if g.transitions[eid].match is not None}
        paths=[]
        for t in terminals:
            frames=frames_for(g,t,traces);mapping=dict(g.states[t].mapping)
            complete=len(mapping)==len(raw['reactant']['elements'])
            events=bond_events(raw,mapping,c.graph_floor,event_tolerance,metal_event_tolerance)
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
        selected_states=set(terminals)|{node for eid in ids
            for node in (g.transitions[eid].source,g.transitions[eid].target)}
        nodes=[s for s in g.states if s.id in selected_states]
        runs.append(dict(label=selection.get('label',f'Context {context_id} · R{c.seed_order[0]} first'),
            context=context_id,seed_order=list(c.seed_order),cuts=list(c.cuts),paths=paths,
            default_terminal=selection.get('default_terminal',terminals[0]),focus=selection.get('focus'),
            config=dict(graph_floor=c.graph_floor,iso_tolerance=c.iso_tolerance,branch_limit=c.branch_limit,
                        self_tolerance=c.self_tolerance,
                        growth_policy=getattr(c,"growth_policy","legacy_wbo_priority"),
                        random_seed=getattr(c,"random_seed",42),
                        event_tolerance=event_tolerance,metal_event_tolerance=metal_event_tolerance),
            graph=dict(states=[dict(id=s.id,assigned=len(s.mapping)) for s in nodes],
                transitions=[dict(id=eid,source=g.transitions[eid].source,target=g.transitions[eid].target,
                    seed=g.transitions[eid].seed,size=len(g.transitions[eid].match.source_atoms) if g.transitions[eid].match else 0)
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
              'Growth frames show one feasible witness; saturation shows canonical placements with '
              'their local relations and shared global placement choices. '
              'Every fragment and placement relationship is checked against its archive. '
              'The endpoint coordinates stay fixed throughout playback.',checks=checks)


def with_represented_targets(document):
    """Project exact per-atom choices from local and shared placement relations.

    The viewer displays marginals only. They are never independent atom choices
    and are computed from the same complete relation used by decoding.
    """
    from copy import deepcopy
    import z3
    from mappa.mappa_core.core.matching.representation.records import Fragment
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation

    doc, cache = deepcopy(document), {}
    for run in doc['runs']:
        for path in run['paths']:
            for frame in path['frames']:
                generators = tuple(tuple(g) for g in frame.get('placement_generators', ()))
                for record in frame['candidates']:
                    if 'relation' not in record:
                        record['target_options'] = {str(r): [p] for r, p in record['witness'].items()}
                        continue
                    fragment = Fragment.from_record(record['relation'])
                    key = fragment.key, generators
                    if key not in cache:
                        system = FiniteConstraints()
                        compiled = compile_relaxed_relation(system, fragment.relation, 'viewer')
                        images = dict(zip(compiled.images, system.act(
                            list(compiled.images.values()), generators, 'viewer_placement')))
                        system.propagate()
                        solver = z3.Solver()
                        system.emit(solver)
                        options = {}
                        for atom, (value, _) in images.items():
                            allowed = []
                            for target in sorted(system.domain(value)):
                                solver.push()
                                solver.add(system.predicate(value, target))
                                status = solver.check()
                                if status == z3.unknown:
                                    raise RuntimeError('fragment target projection is unknown')
                                if status == z3.sat:
                                    allowed.append(target)
                                solver.pop()
                            options[str(atom)] = allowed
                        cache[key] = options
                    record['target_options'] = cache[key]
    return doc
