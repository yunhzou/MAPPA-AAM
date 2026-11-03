"""Existence-only growth must recover the same compressed result at completion."""
from itertools import permutations
from unittest.mock import patch
import os

import numpy as np
import pytest

from mappa import FragmentMatchConfig, match_fragment
from mappa.mappa_core.core.matching.representation.records import FragmentMatchContext
from mappa.mappa_core.core.matching.policy import BondTolerance, DEFAULT_NODE_POLICY
from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from growth_cases import graph as make_graph, represented, discovery_mappings, placed, group_actions, rings, ring_mappings


def graph(elements, edges):
    return make_graph(len(elements), edges, elements)


def record(result):
    return (result.capped, result.branch_count, result.branch_limit,
            [fragment.to_record() for fragment in result.matches], result.deferred_edges,
            set(group_actions(result.placement_generators, len(result.placement_generators[0])))
            if result.placement_generators else set())


def explicit(source, target, atoms, tolerance, edges):
    """Direct tiny-case truth table, independent of production matching."""
    result = set()
    for images in permutations(target, len(atoms)):
        mapping = dict(zip(atoms, images))
        if any(source.nodes[r]['element'] != target.nodes[p]['element']
               for r, p in mapping.items()):
            continue
        if all(target.has_edge(mapping[a], mapping[b]) and
               abs(source.graph['wbo_matrix'][a,b] -
                   target.graph['wbo_matrix'][mapping[a],mapping[b]]) <= tolerance
               for a, b in edges):
            result.add(tuple(sorted(mapping.items())))
    return result


def compare(source, target, *, context=None, config=None):
    config = config or FragmentMatchConfig(iso_tolerance=.5)
    context = context or FragmentMatchContext()
    profile = []
    after = match_fragment(source, target, seed=0, context=context,
                           config=config, profile=profile)
    other = '0' if os.environ.get('MAPPA_NATIVE') == '1' else '1'
    if config.node_policy is None:
        with patch.dict(os.environ, MAPPA_NATIVE=other):
            reference = match_fragment(source, target, seed=0, context=context, config=config)
        assert record(after) == record(reference)
    if profile and after.matches:
        assert profile[0]['max_cands_before'] == 1
        assert profile[0]['max_cands_after'] == 1
    return after, profile


@pytest.mark.parametrize('native',['0','1'])
@pytest.mark.parametrize('size',[3,12])
@pytest.mark.smoke
def test_ring_growth_preserves_all_placements(native,size,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    ring=graph(['C']*size,[(i,(i+1)%size,1) for i in range(size)])
    result,profile=compare(ring,ring)
    expected={tuple((r,(shift+direction*r)%size) for r in range(size))
              for shift in range(size) for direction in (-1,1)}
    assert len(result.matches)==1
    assert represented(result)==expected
    assert profile[0]['feasibility_calls']==size-1
    assert profile[0]['feasibility_searches'] == 0
    assert profile[0]['completion_generated'] == len(result.matches)


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.parametrize('different_environments', [False, True])
def test_completion_returns_each_ring_atom_set_once(native, different_environments, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    ring_edges=[(i,(i+1)%6,1.) for i in range(6)]
    source=graph(['C']*6,ring_edges)
    target_edges=ring_edges+[(a+6,b+6,w) for a,b,w in ring_edges]
    elements=['C']*12
    if different_environments:
        # Distinct external atoms destroy whole-target exchange symmetry while
        # leaving the two internal ring placements unchanged.
        elements+=['F','Cl']
        target_edges += [(0,12,1.),(6,13,1.)]
    target=graph(elements,target_edges)
    result=match_fragment(source,target,seed=0,
        config=FragmentMatchConfig(branch_limit=2,iso_tolerance=.5))
    assert not result.capped
    assert len(result.matches) == (2 if different_environments else 1)
    assert {tuple(sorted(p for _, p in m)) for m in represented(result)} == {
        tuple(range(6)), tuple(range(6,12))}
    for fragment in result.matches:
        assert tuple(fragment.source_atoms)==tuple(range(6))
        assert len(fragment.accepted_edges)==5
    if different_environments:
        assert result.matches[0].source_self_automorph is result.matches[1].source_self_automorph
    capped=match_fragment(source,target,seed=0,
        config=FragmentMatchConfig(branch_limit=1,iso_tolerance=.5))
    assert capped.capped == different_environments
    assert capped.branch_count == (2 if different_environments else 1)
    assert bool(capped.matches) != different_environments


@pytest.mark.parametrize('native', ['0', '1'])
def test_completion_preserves_overlapping_images_and_locked_target_occupancy(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source=graph(['C']*4,[(0,1,1.),(1,2,1.)])
    target=graph(['C']*4,[(0,1,1.),(1,2,1.),(2,3,1.)])
    config=FragmentMatchConfig(branch_limit=2,iso_tolerance=.5)
    result=match_fragment(source,target,seed=0,config=config)
    assert len(result.matches) == 1
    assert {tuple(sorted(p for _, p in m)) for m in represented(result)} == {
        (0,1,2), (1,2,3)}
    for fragment in result.matches:
        assert tuple(fragment.source_atoms)==(0,1,2)
        assert tuple(fragment.accepted_edges)==((0,1),(1,2))
    locked=match_fragment(source,target,seed=0,config=config,
        context=FragmentMatchContext(locked_mapping={3:0},islands={3:1}))
    assert len(locked.matches)==1
    assert tuple(locked.matches[0].target_atoms)==(1,2,3)


@pytest.mark.parametrize('native',['0','1'])
def test_ring_hydrogens_stay_compressed(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    source,target=rings()
    result,_=compare(source,target,config=FragmentMatchConfig(iso_tolerance=1))
    assert len(result.matches) == 1
    assert len({frozenset(p for _, p in m) for m in represented(result)}) == 64
    assert represented(result)==ring_mappings()


@pytest.mark.parametrize('native',['0','1'])
def test_locks_tolerance_atom_sets_components_and_caps(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    source=graph(['X','C','C'],[(0,1,1),(0,2,1)])
    target=graph(['X','C','C','C'],[(0,1,1),(0,2,1.1),(0,3,1.2)])
    result,_=compare(source,target,config=FragmentMatchConfig(branch_limit=3))
    assert not result.capped and len(result.matches)==3
    assert len(represented(result))==6
    source=graph(['O','H','C'],[(0,1,1),(0,2,1)])
    target=graph(['O','H','C']*2,[(0,1,1),(0,2,1),(3,4,1),(3,5,1)])
    result,_=compare(source,target,config=FragmentMatchConfig(branch_limit=1))
    assert not result.capped and result.branch_count == 1
    assert len(represented(result)) == 2
    compare(source,target,config=FragmentMatchConfig(branch_limit=2))
    compare(graph(['X','C','O','N'],[(0,1,1),(1,2,1),(2,3,1),(1,3,1)]),
            graph(['X','C','C','O','N'],[(0,1,1),(0,2,1),(1,3,1),(2,3,1),(3,4,1),(2,4,1)]),
            context=FragmentMatchContext(locked_mapping={2:3,3:4},islands={2:1,3:1}))
    compare(graph(['C','O','O'],[(0,1,1)]),graph(['C','O','O'],[(0,1,1)]),
            context=FragmentMatchContext(locked_mapping={2:1},islands={2:1}))
    anchor=graph(['C','O'],[(0,1,1)])
    compare(anchor,anchor,context=FragmentMatchContext(locked_mapping={0:0},islands={0:1}),
            config=FragmentMatchConfig(allow_mapped_seed=True))
    for inclusive in (True,False):
        compare(graph(['C','O'],[(0,1,1)]),graph(['C','O'],[(0,1,1.25)]),
                config=FragmentMatchConfig(bond_tolerance=BondTolerance(.25,inclusive=inclusive)))


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.parametrize('limit', [1, 2, 3])
def test_completed_choices_obey_the_branch_limit(native, limit, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['O', 'H', 'C'], [(0, 1, 1), (0, 2, 1)])
    target = graph(['O', 'H', 'C'] * 2 + ['F', 'Cl'],
                   [(0, 1, 1), (0, 2, 1), (3, 4, 1), (3, 5, 1),
                    (2, 6, 1), (5, 7, 1)])
    result, _ = compare(source, target, config=FragmentMatchConfig(branch_limit=limit))
    assert result.branch_count == 2
    assert result.capped == (limit < 2)


@pytest.mark.parametrize('native',['0','1'])
def test_growth_preserves_only_accepted_edges_and_can_reassign(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    ring=graph(['C']*3,[(0,1,1),(1,2,1),(0,2,1)])
    chain=graph(['C']*3,[(0,1,1),(1,2,1)])
    result,_=compare(ring,chain)
    assert len(result.matches[0].source_atoms)==3 and len(result.matches[0].accepted_edges)==2
    compare(chain,ring)
    _, profile = compare(graph(['C','C','O'],[(0,1,1),(1,2,1)]),
                         graph(['C','C','O'],[(0,1,1),(0,2,1)]))
    assert profile[0]['feasibility_searches'] > 0
    result,_=compare(graph(['C','O','N'],[(0,1,1),(0,2,1)]),
            graph(['C','N'],[(0,1,1)]),config=FragmentMatchConfig(random_seed=2))
    assert result.matches[0].source_atoms==(0,2)
    assert result.deferred_edges==frozenset({(0,1)})


def test_query_can_change_a_preferred_assignment():
    source=graph(['C','C','O'],[(0,1,1),(1,2,1)])
    target=graph(['C','C','O'],[(0,1,1),(0,2,1)])
    query=FragmentQuery(source,target,{0,1,2},[(0,1),(1,2)],{},DEFAULT_NODE_POLICY,.2)
    assert query.find(preferred={0:0,1:1})=={0:1,1:0,2:2}


@pytest.mark.parametrize('native',['0','1'])
def test_completion_uses_constraints_beyond_the_prefix(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    source=graph(['C','C','O'],[(0,1,1),(1,2,1)])
    target=graph(['C']*6+['O'],[(i,(i+1)%6,1) for i in range(6)]+[(0,6,1)])
    result,profile=compare(source,target)
    assert len(result.matches) == profile[0]['completion_generated'] == 1
    assert represented(result) == {
        ((0, neighbor), (1, 0), (2, 6)) for neighbor in (1, 5)}


@pytest.mark.parametrize('native',['0','1'])
def test_random_weighted_graphs_recover_complete_records(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    rng=np.random.default_rng(1031)
    for case in range(80):
        graphs=[]
        for side in range(2):
            size=4+(side==1 and case%3==0)
            elements=['C']*size if case%2==0 else list(rng.choice(['C','O','N'],size))
            edges=[(a,b,float(w)) for a in range(size) for b in range(a+1,size)
                   if (w:=rng.choice([0.,.1,.8,1.,1.2,1.5]))]
            graphs.append(graph(elements,edges))
        result, _ = compare(*graphs,config=FragmentMatchConfig(random_seed=case,branch_limit=10000))
        if result.branch_count > 1:
            for limit in (result.branch_count - 1, result.branch_count):
                compare(*graphs, config=FragmentMatchConfig(random_seed=case, branch_limit=limit))


def test_trace_identifies_witnesses_and_completion(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE','0')
    ring=graph(['C']*3,[(0,1,1),(1,2,1),(0,2,1)])
    events=[]
    result=match_fragment(ring,ring,seed=0,
        config=FragmentMatchConfig(),events=events)
    start=next(e for e in events if e['type']=='seed_start')
    assert start['candidate_kind']=='witness'
    assert all(e['cands_after']==1 for e in events if e['type']=='commit')
    assert len(result.matches)==1


def test_removed_growth_options_are_rejected():
    from mappa import AAMSearchConfig
    from mappa.mappa_core.core.matching.growth import native
    assert not hasattr(native._engine, 'GrowthReplay')
    for config in (FragmentMatchConfig, AAMSearchConfig):
        with pytest.raises(TypeError, match='growth_policy'):
            config(growth_policy='incremental')


def test_fragment_cache_requires_a_completed_result_view():
    ring = graph(['C'] * 3, [(0, 1, 1), (1, 2, 1), (0, 2, 1)])
    with pytest.raises(TypeError, match='fragment cache'):
        match_fragment(ring, ring, seed=0,
            context=FragmentMatchContext(fragment_cache=object()))


@pytest.mark.parametrize('native', ['0', '1'])
def test_distinct_per_bond_tolerances(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source = graph(['O', 'C', 'C', 'C'], [(0, i, 1.) for i in range(1, 4)])
    target = graph(['O', 'C', 'C', 'C'], [(0, 1, .875), (0, 2, 1.125), (0, 3, 1.375)])
    limits = {1: .2, 2: .4, 3: .2}
    tolerance = BondTolerance(lambda g, a, b: limits[max(a, b)], inclusive=False)
    result, _ = compare(source, target, config=FragmentMatchConfig(bond_tolerance=tolerance))
    assert len(set().union(*(discovery_mappings(fragment) for fragment in result.matches))) == 2
    assert len(represented(result)) == 6


def test_custom_node_attributes(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(['C'] * 3, [(0, 1, 1), (1, 2, 1)])
    target = graph(['C'] * 3, [(0, 1, 1), (0, 2, 1)])
    for g in (source, target):
        for atom in g:
            g.nodes[atom]['label'] = 'end' if atom == 2 else 'any'
    result, _ = compare(source, target, config=FragmentMatchConfig(node_policy=('element', 'label')))
    assert dict(result.matches[0].witness) == {0: 1, 1: 0, 2: 2}


def test_all_alternatives_are_constructed_only_after_growth(monkeypatch):
    from mappa.mappa_core.core.matching.growth import engine
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    ring = graph(['C'] * 6, [(i, (i + 1) % 6, 1) for i in range(6)])
    events = []
    construct = engine.complete_fragment

    def complete(*args, **kwargs):
        assert sum(e['type'] == 'commit' for e in events) == 5
        return construct(*args, **kwargs)

    monkeypatch.setattr(engine, 'complete_fragment', complete)
    result = match_fragment(ring, ring, seed=0, events=events,
        config=FragmentMatchConfig())
    assert len(represented(result)) == 12


@pytest.mark.parametrize('native', ['0', '1'])
def test_random_locked_queries_preserve_the_complete_relation(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    rng = np.random.default_rng(81031)
    for case in range(40):
        graphs = []
        for size in (4, 5):
            edges = [(a, b, float(w)) for a in range(size) for b in range(a + 1, size)
                     if (w := rng.choice([0., .1, .8, 1., 1.2, 1.5]))]
            graphs.append(graph(['C'] * size, edges))
        source, target = graphs
        context = FragmentMatchContext(locked_mapping={3: 4}, islands={3: 1})
        config = FragmentMatchConfig(random_seed=case, branch_limit=10000)
        result, _ = compare(source, target, context=context, config=config)
        if result.matches:
            fragment = result.matches[0]
            expected = explicit(source, target, fragment.source_atoms, .5, fragment.accepted_edges)
            expected = {m for m in expected if all(p != 4 or r == 3 for r, p in m)
                        and (3 not in dict(m) or dict(m)[3] == 4)}
            discovered = placed(set().union(*(discovery_mappings(match) for match in result.matches)),
                                result.placement_generators)
            assert discovered == expected
            expanded = represented(result)
            assert discovered <= expanded
            assert all(p != 4 or r == 3 for mapping in expanded for r, p in mapping)


def test_native_witness_observation_preserves_a_locked_seed():
    from mappa.mappa_core.core.matching.growth import native
    source = graph(['C', 'O'], [(0, 1, 1.)])
    source_view = native.source_graph(source)
    target_view = native.target_graph(source)
    result = native._engine.grow_island(source_view.graph, target_view.graph,
        0, [0, -1], .2, .5, 1, 100, [(0, 1)], [], True,
        observe=True)
    assert result['observation']['candidate_count'] == 1
    assert result['observation']['witness'][0] == 0
    assert result['observation']['anchors'] == [(0, 0)]


def normalized(value):
    if isinstance(value, dict):
        return {str(k): normalized(v) for k, v in value.items()}
    if isinstance(value, (set, frozenset)):
        return [normalized(v) for v in sorted(value)]
    if isinstance(value, (tuple, list)):
        return [normalized(v) for v in value]
    return value


@pytest.mark.parametrize('native', ['0', '1'])
def test_saved_growth_constraints_and_mapping_coverage(native, monkeypatch):
    import json
    from pathlib import Path
    monkeypatch.setenv('MAPPA_NATIVE', native)
    fixture = json.loads((Path(__file__).parent / 'data/accepted_growth_records.json').read_text())
    for case in fixture['cases']:
        source, target = [graph(**value) for value in case['inputs']]
        result = match_fragment(source, target, seed=0, config=FragmentMatchConfig(
            random_seed=case['random_seed'], branch_limit=10000))
        expected = case['expected'][3]
        assert bool(result.matches) == bool(expected)
        if not expected:
            continue
        fragment, bonds, deferred = expected[0][1:4]
        assert normalized(result.deferred_edges) == deferred
        for match in result.matches:
            assert list(match.source_atoms) == fragment
            assert normalized(match.accepted_edges) == bonds
        allowed = explicit(source, target, fragment, .5, bonds)
        assert placed(set().union(*(discovery_mappings(match) for match in result.matches)),
                      result.placement_generators) == allowed


def test_completed_relation_propagates_fixed_image_domains_before_saving(monkeypatch):
    """The forced star center cannot leak into any leaf's saved response domain."""
    from growth_cases import graph
    from mappa import match_fragment, FragmentMatchConfig

    star = graph(13, [(0, leaf) for leaf in range(1, 13)])
    saved = []
    for backend in ('0', '1'):
        monkeypatch.setenv('MAPPA_NATIVE', backend)
        fragment, = match_fragment(star, star, seed=1,
            config=FragmentMatchConfig(iso_tolerance=.1)).matches
        domains = dict(fragment.relation.domains)
        assert domains[0] == (0,)
        assert all(domains[leaf] == tuple(range(1, 13)) for leaf in range(1, 13))
        assert fragment.relation.allows_discovery(dict(fragment.witness))
        saved.append(fragment)
    assert saved[0] == saved[1]
