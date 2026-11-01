"""Existence-only growth must recover the same compressed result at completion."""
from itertools import permutations, product
from unittest.mock import patch
import os

from sympy.combinatorics import Permutation, PermutationGroup

import numpy as np
import pytest

from mappa import FragmentMatchConfig, match_fragment
from mappa.mappa_core.core.matching.representation.records import FragmentMatchContext
from mappa.mappa_core.core.matching.representation.state import candidate_from_record
from mappa.mappa_core.core.matching.policy import BondTolerance, DEFAULT_NODE_POLICY
from mappa.mappa_core.core.matching.symmetry.canonical import _CandidateAutomorphismCanonicalizer
from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from test_growth_compression_exactness import graph, represented
from test_target_action_growth import rings, ring_mappings


def record(result):
    return (result.capped, result.branch_count, result.branch_limit,
            [(dict(m), m.fragment, m.preserved_bonds, m.deferred_edges, m.symmetry)
             for m in result.matches])


def relation(result, target, config, context):
    """Independent small-case expansion of pools and the exact conditioned group."""
    canonicalizer = _CandidateAutomorphismCanonicalizer(target,
        wbo_tol=config.iso_tolerance, locked_mapping=context.locked_mapping,
        node_policy=config.node_policy)
    represented = set()
    for match in result.matches:
        candidate = candidate_from_record(match.symmetry)
        generators = canonicalizer.atom_generators(candidate)
        group = PermutationGroup([Permutation(list(g)) for g in generators]
            or [Permutation(list(range(len(target))))])
        choices = [permutations(block.p_atoms, len(block.r_atoms)) for block in candidate.blocks]
        for selection in product(*choices):
            mapping = dict(candidate.mapping)
            for block, images in zip(candidate.blocks, selection):
                mapping.update(zip(block.r_atoms, images))
            for action in group.generate_schreier_sims():
                represented.add(tuple(sorted((r, action(p)) for r, p in mapping.items())))
    return represented


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
    candidates=[candidate_from_record(m.symmetry) for m in result.matches]
    expected={tuple((r,(shift+direction*r)%size) for r in range(size))
              for shift in range(size) for direction in (-1,1)}
    assert len(result.matches)==1
    assert represented(candidates,ring,.5)==expected
    assert profile[0]['feasibility_calls']==size-1
    assert profile[0]['feasibility_searches'] == 0
    assert profile[0]['completion_generated'] == len(result.matches)


@pytest.mark.parametrize('native',['0','1'])
def test_ring_hydrogens_stay_compressed(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    source,target=rings()
    result,_=compare(source,target,config=FragmentMatchConfig(iso_tolerance=1))
    assert len(result.matches)==1
    assert represented([candidate_from_record(m.symmetry) for m in result.matches],target,1)==ring_mappings()


@pytest.mark.parametrize('native',['0','1'])
def test_locks_tolerance_pools_components_and_caps(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    source=graph(['X','C','C'],[(0,1,1),(0,2,1)])
    target=graph(['X','C','C','C'],[(0,1,1),(0,2,1.1),(0,3,1.2)])
    result,_=compare(source,target,config=FragmentMatchConfig(branch_limit=1))
    assert not result.capped and len(result.matches)==1
    assert len(represented([candidate_from_record(result.matches[0].symmetry)],target,.5))==6
    source=graph(['O','H','C'],[(0,1,1),(0,2,1)])
    target=graph(['O','H','C']*2,[(0,1,1),(0,2,1),(3,4,1),(3,5,1)])
    result,_=compare(source,target,config=FragmentMatchConfig(branch_limit=1))
    assert result.capped and result.branch_count==2
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
    target = graph(['O', 'H', 'C'] * 2,
                   [(0, 1, 1), (0, 2, 1), (3, 4, 1), (3, 5, 1)])
    result, _ = compare(source, target, config=FragmentMatchConfig(branch_limit=limit))
    assert result.branch_count == 2
    assert result.capped == (limit < 2)


@pytest.mark.parametrize('native',['0','1'])
def test_growth_preserves_only_accepted_edges_and_can_reassign(native,monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE',native)
    ring=graph(['C']*3,[(0,1,1),(1,2,1),(0,2,1)])
    chain=graph(['C']*3,[(0,1,1),(1,2,1)])
    result,_=compare(ring,chain)
    assert len(result.matches[0])==3 and len(result.matches[0].preserved_bonds)==2
    compare(chain,ring)
    _, profile = compare(graph(['C','C','O'],[(0,1,1),(1,2,1)]),
                         graph(['C','C','O'],[(0,1,1),(0,2,1)]))
    assert profile[0]['feasibility_searches'] > 0
    result,_=compare(graph(['C','O','N'],[(0,1,1),(0,2,1)]),
            graph(['C','N'],[(0,1,1)]),config=FragmentMatchConfig(random_seed=2))
    assert result.matches[0].fragment==frozenset({0,2})
    assert result.matches[0].deferred_edges==frozenset({(0,1)})


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
    assert represented([candidate_from_record(result.matches[0].symmetry)], target, .5) == {
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
    assert len(represented([candidate_from_record(m.symmetry) for m in result.matches], target, .5)) == 2


def test_custom_node_attributes(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(['C'] * 3, [(0, 1, 1), (1, 2, 1)])
    target = graph(['C'] * 3, [(0, 1, 1), (0, 2, 1)])
    for g in (source, target):
        for atom in g:
            g.nodes[atom]['label'] = 'end' if atom == 2 else 'any'
    result, _ = compare(source, target, config=FragmentMatchConfig(node_policy=('element', 'label')))
    assert result.matches[0].mapping == {0: 1, 1: 0, 2: 2}


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
    assert len(represented([candidate_from_record(m.symmetry) for m in result.matches], ring, .5)) == 12


@pytest.mark.parametrize('native', ['0', '1'])
def test_random_locked_queries_preserve_the_complete_relation(native, monkeypatch):
    from test_growth_compression_exactness import explicit
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
            expected = explicit(source, target, sorted(fragment.fragment), .5, fragment.preserved_bonds)
            expected = {m for m in expected if all(p != 4 or r == 3 for r, p in m)
                        and (3 not in dict(m) or dict(m)[3] == 4)}
            assert relation(result, target, config, context) == expected


def test_native_witness_observation_preserves_a_locked_seed():
    from mappa.mappa_core.core.matching.growth import native
    from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits
    source = graph(['C', 'O'], [(0, 1, 1.)])
    source_view = native.source_graph(source)
    target_view = native.target_graph(source, _nauty_orbits(source, wbo_tol=.5))
    result = native._engine.grow_island(source_view.graph, target_view.graph,
        0, [0, -1], .2, .5, 1, 100, [(0, 1)], [], True,
        observe=True)
    assert result['observation']['candidate_count'] == 1
    assert result['observation']['symmetry']['witness'] == {0: 0}
    assert result['observation']['symmetry']['exact_fixed'] == [0]


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
    from test_growth_compression_exactness import explicit
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
        old = [candidate_from_record(match[4]) for match in expected]
        new = [candidate_from_record(match.symmetry) for match in result.matches]
        fragment, bonds, deferred = expected[0][1:4]
        for match in result.matches:
            assert sorted(match.fragment) == fragment
            assert normalized(match.preserved_bonds) == bonds
            assert normalized(match.deferred_edges) == deferred
        allowed = explicit(source, target, fragment, .5, bonds)
        assert represented(old, target, .5) == represented(new, target, .5) == allowed
