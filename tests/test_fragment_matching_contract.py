"""Release gates for one-witness existence and distinct-atom-set completion."""
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import resource
import sys
import time
import traceback
from unittest.mock import patch

import networkx as nx
import pytest

from growth_cases import graph, represented
from mappa import BondTolerance, FragmentMatchConfig, match_fragment
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from mappa.mappa_core.core.matching.policy import DEFAULT_NODE_POLICY
from mappa.mappa_core.core.matching.representation.records import Fragment


FIXTURE = Path(__file__).parent / 'fixtures/fragment_matching/enzymemap_16175.json'


def rings(copies=1):
    return graph(6 * copies, [(6 * c + i, 6 * c + (i + 1) % 6)
                             for c in range(copies) for i in range(6)])


def saved_case():
    saved = json.loads(FIXTURE.read_text())
    endpoints = [graph(len(saved[side]['elements']), saved[side]['weights'],
                       saved[side]['elements']) for side in ('source', 'target')]
    return *endpoints, saved


def query_case(name):
    if name in ('slow_yes', 'real_no'):
        source, target, saved = saved_case()
        return source, target, {**saved['queries'][name], 'tolerance': saved['tolerance']}
    if name == 'two_rings':
        source, target = rings(), rings(2)
        return source, target, dict(atoms=list(source), bonds=list(source.edges),
                                   locked=[], tolerance=.2)
    if name in ('articulation_yes', 'articulation_no'):
        count = 13 if name.endswith('_yes') else 14
        source = graph(count, [(r, r + 1) for r in range(count - 1)])
        target = graph(19, [(0, 1 + 6 * c) for c in range(3)]
            + [(1 + 6 * c + i, 1 + 6 * c + (i + 1) % 6)
               for c in range(3) for i in range(6)])
        return source, target, dict(atoms=list(source), bonds=list(source.edges),
                                   locked=[], tolerance=.2)
    # All carbon atoms pass element/degree tests and the total target pool is
    # large enough. Twelve of them compete for only eleven reachable images.
    count = 12 if name == 'competing_no' else 11
    source = graph(count + 3, [(0, r) for r in range(2, count + 2)] + [(1, count + 2)],
                   ['O', 'N'] + ['C'] * (count + 1))
    target = graph(34, [(0, p) for p in range(2, 13)] + [(0, 33)]
                   + [(1, p) for p in range(13, 33)], ['O', 'N'] + ['C'] * 31 + ['H'])
    return source, target, dict(atoms=list(source), bonds=list(source.edges),
                               locked=[(0, 0), (1, 1)], tolerance=.2)


def assert_witness(source, target, query, witness):
    """Check the returned assignment directly, without a matcher or decoder."""
    assert set(witness) == set(query['atoms'])
    assert len(set(witness.values())) == len(witness)
    assert all(source.nodes[r]['element'] == target.nodes[p]['element']
               for r, p in witness.items())
    assert all(witness[r] == p for r, p in query['locked'] if r in witness)
    occupied = {p for r, p in query['locked'] if r not in witness}
    assert not occupied.intersection(witness.values())
    for a, b in query['bonds']:
        p, q = witness[a], witness[b]
        assert target.has_edge(p, q)
        assert abs(source.graph['wbo_matrix'][a, b] - target.graph['wbo_matrix'][p, q]) <= query['tolerance']


def matching_witnesses(source, target, query, bond_tolerance=None):
    """Compare the accelerated kernel with the independent Python solver."""
    locked = dict(query['locked'])
    source_view = native.source_graph(source, bond_tolerance)
    target_view = native.target_graph(target)
    result = native._engine.find_fragment_witness(
        source_view.graph, target_view.graph, query['atoms'], query['bonds'],
        [locked.get(r, -1) for r in source], query['tolerance'])['witness']
    expected = FragmentQuery(source, target, query['atoms'], query['bonds'],
        locked, DEFAULT_NODE_POLICY, query['tolerance'], bond_tolerance).find()
    assert result == expected
    return result


@pytest.mark.parametrize('connected', [False, True], ids=['absent', 'present'])
def test_sparse_locked_query_ignores_unavailable_targets(connected):
    size = 65
    source = graph(4, [(0, 1), (0, 2), (1, 2)], ['O', 'C', 'C', 'C'])
    occupied, first, last = range(size - 3, size)
    edges = [(0, p) for p in (occupied, first, last)] + [(occupied, first), (occupied, last)]
    if connected:
        edges.append((first, last))
    target = graph(size, edges, ['O'] + ['He'] * (size - 4) + ['C'] * 3)
    query = dict(atoms=[0, 1, 2], bonds=list(source.edges),
                 locked=[(0, 0), (3, occupied)], tolerance=.2)
    result = matching_witnesses(source, target, query)
    # Only these two free carbons can fill the triangle. Both pass degree
    # filtering even when their connecting edge is absent; the occupied third
    # carbon cannot rescue that missing edge.
    expected = {((0, 0), (1, a), (2, b)) for a, b in ((first, last), (last, first))}
    if connected:
        assert tuple(sorted(result.items())) in expected
        assert_witness(source, target, query, result)
    else:
        assert result is None


@pytest.mark.parametrize('inclusive', [False, True])
@pytest.mark.parametrize('position', [-1, 0, 1])
def test_conditional_tolerance_support_at_exact_threshold(inclusive, position):
    source = graph(4, [(0, 1), (0, 2), (1, 2)], ['O', 'C', 'C', 'C'])
    weight = 1.25
    if position:
        weight = math.nextafter(weight, -math.inf if position < 0 else math.inf)
    target = graph(65, [(0, 62), (0, 63, weight), (0, 64, 1.125), (63, 64)],
                   ['O'] + ['He'] * 61 + ['C'] * 3)
    limits = {(0, 1): .25, (0, 2): .5, (1, 2): .125}
    policy = BondTolerance(lambda graph, a, b: limits[a, b], inclusive=inclusive).resolve(source)
    query = dict(atoms=[0, 1, 2], bonds=list(source.edges),
                 locked=[(0, 0), (2, 64), (3, 62)], tolerance=.5)
    result = matching_witnesses(source, target, query, policy)
    difference = abs(weight - 1.)
    allowed = difference <= .25 if inclusive else difference < .25
    assert result == ({0: 0, 1: 63, 2: 64} if allowed else None)


def query_worker(engine, name, connection):
    """Run the real existence kernel with a killable external time guard."""
    try:
        source, target, query = query_case(name)
        locked = dict(query['locked'])
        if engine == 'native':
            source_view = native.source_graph(source)
            target_view = native.target_graph(target)
            started = time.perf_counter()
            result = native._engine.find_fragment_witness(
                source_view.graph, target_view.graph, query['atoms'], query['bonds'],
                [locked.get(r, -1) for r in source], query['tolerance'])
        else:
            started = time.perf_counter()
            solver = FragmentQuery(source, target, query['atoms'], query['bonds'],
                                   locked, DEFAULT_NODE_POLICY, query['tolerance'])
            result = dict(witness=solver.find(), nodes=solver.nodes)
        result['seconds'] = time.perf_counter() - started
        result['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
            1 if sys.platform == 'darwin' else 1024)
        connection.send(result)
    except Exception:
        connection.send(dict(error=traceback.format_exc()))
    finally:
        connection.close()


def bounded_query(engine, name, worker=query_worker):
    context = mp.get_context('spawn')
    receiver, sender = context.Pipe(duplex=False)
    process = context.Process(target=worker, args=(engine, name, sender))
    process.start()
    sender.close()
    try:
        assert receiver.poll(10), 'existence query exceeded 10 s; this is a failure, not a no-match result'
        result = receiver.recv()
        assert 'error' not in result, result.get('error')
        return result
    finally:
        process.join(1)
        if process.is_alive():
            process.terminate()
            process.join()
        receiver.close()


@pytest.mark.parametrize('engine', ['python', 'native'])
@pytest.mark.parametrize('name', ['two_rings', 'slow_yes', 'real_no', 'competing_yes', 'competing_no'])
def test_same_existence_function_finds_one_or_proves_none(engine, name, record_property):
    source, target, query = query_case(name)
    result = bounded_query(engine, name)
    record_property('measurement', json.dumps({**result, 'witness': result['witness'] is not None}))
    print(engine, name, {**result, 'witness': result['witness'] is not None})
    if name.endswith('_no'):
        assert result['witness'] is None
        assert result['nodes'] == 1, 'this contradiction must be proved before branching on assignments'
    else:
        assert result['witness'] is not None
        assert_witness(source, target, query, result['witness'])
        assert result['nodes'] <= 2 * len(query['atoms']), 'existence is exploring redundant assignments'
    assert result['seconds'] < 2, 'small fixed regression query exceeded its kernel budget'
    assert result['peak_rss_bytes'] < 512 * 2**20


def test_real_no_match_has_an_independent_proof_and_a_matching_predecessor():
    source, target, saved = saved_case()
    positive, negative = saved['queries']['slow_yes'], saved['queries']['real_no']
    assert_witness(source, target, {**positive, 'tolerance': 1.}, dict(positive['witness']))
    previous = dict(negative['previous_witness'])
    proposal = negative['proposal']
    assert set(negative['atoms']) == set(previous) | {proposal[1]}
    before = dict(atoms=sorted(previous), bonds=[edge for edge in negative['bonds']
                  if set(edge) != set(proposal)], locked=[], tolerance=1.)
    assert_witness(source, target, before, previous)
    # This 28-atom accepted subtree cannot fit even when all bond weights are
    # ignored. NetworkX supplies an independent, stronger absence proof.
    obstruction = negative['obstruction']
    assert {frozenset(e) for e in obstruction['edges']} <= {frozenset(e) for e in negative['bonds']}
    pattern = nx.Graph()
    pattern.add_nodes_from((r, source.nodes[r]) for r in obstruction['atoms'])
    pattern.add_edges_from(obstruction['edges'])
    oracle = nx.isomorphism.GraphMatcher(target, pattern,
        node_match=lambda a, b: a['element'] == b['element'])
    assert not oracle.subgraph_is_monomorphic()


@pytest.mark.parametrize('engine', ['python', 'native'])
def test_articulation_obstruction_has_a_matching_predecessor(engine):
    source, target, query = query_case('articulation_yes')
    witness = dict(enumerate([2, 3, 4, 5, 6, 1, 0, 7, 8, 9, 10, 11, 12]))
    assert_witness(source, target, query, witness)
    without_hub = target.copy()
    without_hub.remove_node(0)
    assert sorted(map(len, nx.connected_components(without_hub))) == [6, 6, 6]
    # A simple path can enter at most two components through a single hub:
    # at most six atoms, then the hub, then six atoms. Fourteen is impossible.
    for name in ('articulation_yes', 'articulation_no'):
        result = bounded_query(engine, name)
        if name.endswith('_yes'):
            assert_witness(source, target, query, result['witness'])
        else:
            assert result['witness'] is None
        assert result['seconds'] < 2
        assert result['peak_rss_bytes'] < 512 * 2**20


def complete_once(engine, source, target, config):
    """Observe completed records directly, independently of profile counters."""
    constructed = []
    profile = []
    initialize = Fragment.__init__

    def observe(fragment, *args, **kwargs):
        initialize(fragment, *args, **kwargs)
        constructed.append((fragment.source_atoms, fragment.target_atoms))

    with patch.dict(os.environ, MAPPA_NATIVE=engine), patch.object(Fragment, '__init__', observe):
        result = match_fragment(source, target, seed=0, config=config, profile=profile)
    assert not result.capped
    assert result.branch_count == len(result.matches)
    assert profile[0]['completion_generated'] == len(result.matches)
    assert constructed == [(fragment.source_atoms, fragment.target_atoms)
                           for fragment in result.matches]
    return result


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
@pytest.mark.parametrize('copies', [1, 2])
def test_full_matching_constructs_one_fragment_per_ring(engine, copies):
    source, target = rings(), rings(copies)
    result = complete_once(engine, source, target,
        FragmentMatchConfig(iso_tolerance=.2, branch_limit=copies))
    assert len(result.matches) == copies
    for component, fragment in enumerate(result.matches):
        assert fragment.source_atoms == tuple(range(6))
        assert fragment.target_atoms == tuple(range(6 * component, 6 * component + 6))
        expected = {tuple((r, 6 * component + (offset + direction * r) % 6) for r in range(6))
                    for offset in range(6) for direction in (-1, 1)}
        assert represented([fragment]) == expected


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_free_hydrogen_keeps_every_target_atom_set(engine):
    source = graph(1, [], ['H'])
    target = graph(10, [(p, p + 5) for p in range(5)],
                   ['H'] * 5 + ['C', 'C', 'O', 'O', 'N'])
    # External attachments do not merge different labelled target atom sets.
    result = complete_once(engine, source, target, FragmentMatchConfig(branch_limit=5))
    assert [fragment.target_atoms for fragment in result.matches] == [(p,) for p in range(5)]
    assert [represented([fragment]) for fragment in result.matches] == [
        {((0, p),)} for p in range(5)]


def asymmetric_pool():
    count = 14
    source = graph(count + 1, [(0, r, 1. + r / 1000) for r in range(1, count + 1)],
                   ['O'] + ['C'] * count)
    target = graph(count + 1,
        [(0, p, 1. + p / 10000) for p in range(1, count + 1)]
        + [(a, b, 1. + (a * count + b) / 10000)
           for a in range(1, count + 1) for b in range(a + 1, count + 1)],
        ['O'] + ['C'] * count)
    return source, target


def completion_worker(engine, name, connection):
    try:
        assert name == 'asymmetric_pool'
        source, target = asymmetric_pool()
        started = time.perf_counter()
        result = complete_once(engine, source, target, FragmentMatchConfig(branch_limit=1))
        connection.send(dict(seconds=time.perf_counter() - started,
            fragment_count=len(result.matches),
            source_atoms=result.matches[0].source_atoms,
            target_atoms=result.matches[0].target_atoms,
            domains=result.matches[0].relation.domains,
            witness=result.matches[0].witness,
            accepted_edges=result.matches[0].accepted_edges,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
                1 if sys.platform == 'darwin' else 1024)))
    except Exception:
        connection.send(dict(error=traceback.format_exc()))
    finally:
        connection.close()


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_large_tolerance_relation_needs_no_bijection_enumeration(engine):
    source, target = asymmetric_pool()
    # Exact WBOs differ, but all fourteen carbon assignments satisfy discovery.
    # Inspect the complete domains without expanding their 14! permutations.
    assert len({target.graph['wbo_matrix'][0, p] for p in range(1, 15)}) == 14
    result = bounded_query(engine, 'asymmetric_pool', worker=completion_worker)
    assert result['fragment_count'] == 1
    assert result['source_atoms'] == result['target_atoms'] == tuple(range(15))
    domains = dict(result['domains'])
    assert domains[0] == (0,)
    assert all(domains[r] == tuple(range(1, 15)) for r in range(1, 15))
    assert set(result['accepted_edges']) == set(source.edges)
    assert_witness(source, target, dict(atoms=list(source), bonds=list(source.edges),
                   locked=[], tolerance=.5), dict(result['witness']))
    assert result['seconds'] < 2
    assert result['peak_rss_bytes'] < 512 * 2**20


@pytest.mark.parametrize('policy_kind', [
    'callable_identity', 'attribute_identity', 'element_subclass_identity', 'callable_pool'])
def test_custom_policy_survives_completed_relations_and_records(policy_kind, monkeypatch):
    from mappa.mappa_core.core.matching.policy import ElementNodeMatchPolicy
    from mappa.mappa_core.domain import AAMHierarchy
    from mappa.mappa_core.core.search.branching.construction import find_islands

    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source = graph(2, [(0, 1)])
    size = 4 if policy_kind == 'callable_pool' else 3
    target = graph(size, [(a, b) for a in range(size) for b in range(a + 1, size)])
    if policy_kind == 'attribute_identity':
        for endpoint in (source, target):
            nx.set_node_attributes(endpoint, {atom: atom for atom in endpoint}, 'role')
        policy = ('element', 'role')
    elif policy_kind == 'callable_pool':
        policy = lambda source, r, target, p: (r == 0) == (p == 0)
    elif policy_kind == 'element_subclass_identity':
        class IdentityPolicy(ElementNodeMatchPolicy):
            def compatible(self, source, r, target, p):
                return r == p
        policy = IdentityPolicy()
    else:
        policy = lambda source, r, target, p: r == p
    images = (1, 2, 3) if policy_kind == 'callable_pool' else (1,)
    expected = {((0, 0), (1, image)) for image in images}

    result = match_fragment(source, target, seed=0,
        config=FragmentMatchConfig(node_policy=policy, branch_limit=len(images)))
    assert not result.capped and len(result.matches) == len(images)
    assert represented(result.matches) == expected
    assert all(fragment.target_self_automorph.generators == () for fragment in result.matches)

    search = find_islands(source, target, [0], iso_tol=.5, node_policy=policy,
                          max_branches=len(images))
    paths = list(search.paths())
    assert len(paths) == len(images)
    restored = [AAMHierarchy.from_record(json.loads(json.dumps(path.hierarchy.to_record())))
                for path in paths]
    assert set().union(*(represented(hierarchy.fragments) for hierarchy in restored)) == expected
    monkeypatch.setenv('MAPPA_NATIVE', '1')
    with pytest.raises(ValueError, match='requires the element node policy'):
        match_fragment(source, target, seed=0,
            config=FragmentMatchConfig(node_policy=policy))


def test_custom_policy_actions_preserve_later_fragments(monkeypatch):
    from mappa.mappa_core.core.search.branching.construction import find_islands

    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source, target = graph(2, []), graph(2, [(0, 1)])
    policy = lambda source, r, target, p: r == 0 or p == 1
    search = find_islands(source, target, [0, 1], iso_tol=.5, node_policy=policy)
    complete = [path for path in search.paths() if len(path.mapping) == 2]
    assert len(complete) == 1
    assert complete[0].mapping == {0: 0, 1: 1}
    # Local choices must retain the custom policy after later fragments lock.
    for fragment in complete[0].hierarchy.fragments:
        assert represented([fragment]) == {fragment.witness}
        assert all(policy(source, r, target, p) for r, p in fragment.witness)
