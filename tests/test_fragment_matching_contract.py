"""Release gates for one-witness existence and direct compressed completion.

These tests deliberately fail while completion still replays candidate prefixes
or invokes candidate deduplication. Do not mark those failures xfail.
"""
import json
from contextlib import nullcontext
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
from mappa import FragmentMatchConfig, match_fragment
from mappa.mappa_core.core.matching.growth import native
from mappa.mappa_core.core.matching.growth.existence import FragmentQuery
from mappa.mappa_core.core.matching.policy import DEFAULT_NODE_POLICY
from mappa.mappa_core.core.matching.representation.state import _SymCand, candidate_from_record
from mappa.mappa_core.core.matching.symmetry.orbits import _nauty_orbits


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


def query_worker(engine, name, connection):
    """Run the real existence kernel with a killable external time guard."""
    try:
        source, target, query = query_case(name)
        locked = dict(query['locked'])
        if engine == 'native':
            source_view = native.source_graph(source)
            orbits = _nauty_orbits(target, wbo_tol=query['tolerance'])
            target_view = native.target_graph(target, orbits)
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


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
@pytest.mark.parametrize('copies', [1, 2])
def test_full_matching_constructs_one_candidate_per_ring(engine, copies, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', engine)
    source, target = rings(), rings(copies)
    profile = []
    constructed = []
    initialize = _SymCand.__init__

    def observe(candidate, *args, **kwargs):
        initialize(candidate, *args, **kwargs)
        constructed.append(frozenset(candidate.mapping))

    with monkeypatch.context() as construction:
        if engine == '0':
            construction.setattr(_SymCand, '__init__', observe)
        result = match_fragment(source, target, seed=0,
            config=FragmentMatchConfig(iso_tolerance=.2, branch_limit=copies), profile=profile)
    assert not result.capped
    assert result.branch_count == len(result.matches) == copies
    for component, match in enumerate(result.matches):
        assert set(match.mapping) == set(range(6))
        assert set(match.mapping.values()) == set(range(6 * component, 6 * component + 6))
        expected = {tuple((r, 6 * component + (offset + direction * r) % 6) for r in range(6))
                    for offset in range(6) for direction in (-1, 1)}
        assert represented([candidate_from_record(match.symmetry)], target, .2) == expected
    if engine == '0':
        assert constructed == [frozenset(source)] * copies, (
            'candidate construction replayed prefixes or created extra records: '
            f'{list(map(len, constructed))}')
    assert profile[0]['completion_generated'] == copies, (
        'completion must construct only final compressed candidates; '
        f"constructed {profile[0]['completion_generated']} for {copies} placements")


def test_full_matching_never_invokes_candidate_deduplication(monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', '0')
    # Guard the old imported aliases while they exist. Once removed, this
    # continues to exercise the same public operation without importing them.
    def forbidden(*args, **kwargs):
        pytest.fail('fragment matching invoked the retired candidate deduplication stage')

    for name, module in tuple(sys.modules.items()):
        if name.startswith('mappa.mappa_core.core.matching.') and hasattr(module, '_dedup_sym_cands'):
            monkeypatch.setattr(module, '_dedup_sym_cands', forbidden)
    result = match_fragment(rings(), rings(2), seed=0, config=FragmentMatchConfig(iso_tolerance=.2))
    assert len(result.matches) == 2


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_inconsistent_supplied_orbits_cannot_create_duplicate_families(engine, monkeypatch):
    from mappa import FragmentMatchContext

    monkeypatch.setenv('MAPPA_NATIVE', engine)
    source = graph(2, [(0, 1)])
    target = graph(3, [(0, 1), (1, 2), (2, 0)])
    orbits = _nauty_orbits(target, wbo_tol=.2)
    # This unequal split is not preserved by triangle rotations. Treating it
    # as a class action previously emitted the same six mappings three times.
    orbits.update({0: 0, 1: 0, 2: 1})
    with pytest.raises(ValueError, match='orbit partition'):
        match_fragment(source, target, seed=0,
            context=FragmentMatchContext(target_orbits=orbits),
            config=FragmentMatchConfig(iso_tolerance=.2))


def complete_once(engine, source, target, config):
    """Observe real result construction, independently of the profile counter."""
    constructed = []
    profile = []
    initialize = _SymCand.__init__

    def observe(candidate, *args, **kwargs):
        initialize(candidate, *args, **kwargs)
        constructed.append(frozenset(candidate.mapping))

    observation = patch.object(_SymCand, '__init__', observe) if engine == '0' else nullcontext()
    with patch.dict(os.environ, MAPPA_NATIVE=engine), observation:
        result = match_fragment(source, target, seed=0, config=config, profile=profile)
    assert not result.capped
    assert result.branch_count == len(result.matches)
    assert profile[0]['completion_generated'] == len(result.matches)
    if engine == '0':
        assert constructed == [match.fragment for match in result.matches]
    return result


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_free_hydrogen_pools_span_interchangeable_target_components(engine):
    source = graph(1, [], ['H'])
    target = graph(10, [(p, p + 5) for p in range(5)],
                   ['H'] * 5 + ['C', 'C', 'O', 'O', 'N'])
    # No accepted bond constrains this free hydrogen. Component labels must
    # not split its two interchangeable H-C choices or its two H-O choices.
    result = complete_once(engine, source, target, FragmentMatchConfig(branch_limit=3))
    families = [represented([candidate_from_record(match.symmetry)], target, .5)
                for match in result.matches]
    assert families == [{((0, 0),), ((0, 1),)},
                        {((0, 2),), ((0, 3),)}, {((0, 4),)}]


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_tolerance_equivalent_assignments_remain_one_complete_pool(engine):
    source = graph(4, [(0, 1, 1.2), (0, 2, .8), (0, 3, 1.2), (1, 2, .8)],
                   ['C', 'N', 'C', 'N'])
    target = graph(4, [(0, 1, 1.2), (0, 2, 1.), (1, 3, .8), (2, 3, .8)],
                   ['N', 'C', 'C', 'O'])
    result = complete_once(engine, source, target,
        FragmentMatchConfig(random_seed=11, branch_limit=1))
    assert len(result.matches) == 1
    assert all(match.preserved_bonds == ((0, 1), (1, 2)) for match in result.matches)
    expected = {((0, a), (1, 0), (2, b)) for a, b in ((1, 2), (2, 1))}
    candidate = candidate_from_record(result.matches[0].symmetry)
    assert represented([candidate], target, .5) == expected


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_internal_fragment_atom_can_retain_a_pool(engine):
    source = graph(4, [(0, 1, 1.2), (0, 2, 1.), (0, 3, 1.5),
                       (1, 2, 1.5), (1, 3, .8), (2, 3, .8)], ['O', 'N', 'C', 'O'])
    target = graph(4, [(0, 2, 1.2), (0, 3, 1.5), (1, 2, 1.5),
                       (1, 3, 1.5), (2, 3, .1)], ['O', 'C', 'O', 'O'])
    result = complete_once(engine, source, target,
        FragmentMatchConfig(random_seed=55, branch_limit=1))
    assert len(result.matches) == 1
    match = result.matches[0]
    assert match.preserved_bonds == ((0, 2), (0, 3))
    candidate = candidate_from_record(match.symmetry)
    assert [(block.r_atoms, block.p_atoms) for block in candidate.blocks] == [((0,), (2, 3))]
    assert represented([candidate], target, .5) == {
        ((0, 2), (2, 1), (3, 0)), ((0, 3), (2, 1), (3, 0))}


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_completion_depends_on_constraints_not_growth_order(engine):
    source = graph(4, [(0, 1, .8), (0, 2, 1.), (0, 3, 1.),
                       (1, 2, 1.), (1, 3, 1.5), (2, 3, 1.)])
    target = graph(4, [(0, 1, 1.), (0, 2, 1.5), (0, 3, 1.5),
                       (1, 2, .1), (2, 3, .1)])
    expected = {((0, 0), (1, 1), (2, a), (3, b)) for a, b in ((2, 3), (3, 2))}
    for seed in (8, 9):
        result = complete_once(engine, source, target,
            FragmentMatchConfig(random_seed=seed, branch_limit=1))
        assert len(result.matches) == 1
        match = result.matches[0]
        assert match.preserved_bonds == ((0, 1), (0, 2), (0, 3))
        assert represented([candidate_from_record(match.symmetry)], target, .5) == expected


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
            symmetry=result.matches[0].symmetry,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
                1 if sys.platform == 'darwin' else 1024)))
    except Exception:
        connection.send(dict(error=traceback.format_exc()))
    finally:
        connection.close()


@pytest.mark.parametrize('engine', ['0', '1'], ids=['python', 'native'])
def test_large_tolerance_pool_needs_no_bijection_enumeration(engine):
    source, target = asymmetric_pool()
    # The unique oxygen and its fourteen different exact bond weights fix every
    # carbon under global automorphisms. Compression must therefore come from
    # the valid tolerance pool, not a large exact automorphism group.
    assert len({target.graph['wbo_matrix'][0, p] for p in range(1, 15)}) == 14
    result = bounded_query(engine, 'asymmetric_pool', worker=completion_worker)
    candidate = candidate_from_record(result['symmetry'])
    assert candidate.mapping[0] == 0
    assert len(candidate.blocks) == 1
    block = candidate.blocks[0]
    assert block.r_atoms == block.p_atoms == tuple(range(1, 15))
    # This one complete pool represents all 14! injections without expanding
    # them in either the test or the implementation.
    assert_witness(source, target, dict(atoms=list(source), bonds=list(source.edges),
                   locked=[], tolerance=.5), candidate.mapping)
    assert result['seconds'] < 2
    assert result['peak_rss_bytes'] < 512 * 2**20


@pytest.mark.parametrize('policy_kind', [
    'callable_identity', 'attribute_identity', 'element_subclass_identity', 'callable_pool'])
def test_custom_policy_survives_symmetry_finalization_and_records(policy_kind, monkeypatch):
    from mappa.mappa_core.core.matching.policy import ElementNodeMatchPolicy
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry

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
        config=FragmentMatchConfig(node_policy=policy, branch_limit=1))
    assert not result.capped and len(result.matches) == 1
    candidate = candidate_from_record(result.matches[0].symmetry)
    assert candidate.target_generators is not None
    if policy_kind != 'callable_pool':
        assert candidate.target_generators == ()
    assert represented([candidate], target, .5) == expected

    search = find_islands(source, target, [0], iso_tol=.5, node_policy=policy,
                          max_branches=1)
    finalized, _ = finalize_graph_symmetry(search, target, iso_tolerance=.5)
    paths = list(finalized.paths())
    assert len(paths) == 1
    # Typed/JSON export intentionally omits the finalization marker. The saved
    # exact action must still preserve the policy when the record is reopened.
    record = json.loads(json.dumps(paths[0].hierarchy.to_record()))
    restored = candidate_from_record(record['fragments'][0]['symmetry'])
    assert restored.target_generators == candidate.target_generators
    assert represented([restored], target, .5) == expected
    monkeypatch.setenv('MAPPA_NATIVE', '1')
    with pytest.raises(ValueError, match='requires the element node policy'):
        match_fragment(source, target, seed=0,
            config=FragmentMatchConfig(node_policy=policy))


@pytest.mark.parametrize('workspace_kind', ['python', 'native'])
def test_saved_group_is_stabilized_exactly_under_additional_locks(workspace_kind):
    from mappa.mappa_core.core.matching.symmetry.workspace import (
        ConditionedSymmetryWorkspace, SymmetryWorkspace)

    target = graph(3, [(0, 1), (1, 2), (2, 0)])
    # Neither generator fixes atom 1, but their group contains the exchange
    # of 0 and 2. Filtering generators individually would lose this action.
    generators = ((1, 0, 2), (0, 2, 1))
    candidate = _SymCand({0: 0}, target_generators=generators)
    workspace_type = SymmetryWorkspace if workspace_kind == 'python' else ConditionedSymmetryWorkspace
    workspace = workspace_type(target, .5)
    assert workspace.conditioned_generators(candidate, ((9, 1),)) == ((2, 1, 0),)
    assert candidate.target_generators == generators
    assert workspace.conditioned_generators(candidate, ((8, 0), (9, 1))) == ()


def test_custom_policy_actions_preserve_later_fragments(monkeypatch):
    from mappa.mappa_core.core.search.branching.construction import find_islands
    from mappa.mappa_core.core.search.branching.symmetry import finalize_graph_symmetry

    monkeypatch.setenv('MAPPA_NATIVE', '0')
    source, target = graph(2, []), graph(2, [(0, 1)])
    policy = lambda source, r, target, p: r == 0 or p == 1
    search = find_islands(source, target, [0, 1], iso_tol=.5, node_policy=policy)
    finalized, _ = finalize_graph_symmetry(search, target, iso_tolerance=.5)
    complete = [path for path in finalized.paths() if len(path.mapping) == 2]
    assert len(complete) == 1
    assert complete[0].mapping == {0: 0, 1: 1}
    # An earlier fragment's actions also move later fragments. Swapping these
    # targets would make source atom 1 violate its policy, so neither an open
    # pool nor a global action may introduce that swap.
    for fragment in complete[0].hierarchy.fragments:
        assert all(len(domain.p_atoms) <= 1 for domain in fragment.symmetry_domains)
        assert all(action.images == (0, 1) for action in fragment.target_generators)
