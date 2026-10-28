from dataclasses import replace

import networkx as nx
import numpy as np
import pytest

from mappa import AAMSearchConfig
from mappa.adaptive_search import AdaptiveSearchCondition
from mappa.alignment.sweep import cut_sweep_items
from mappa.cut_replay import FragmentRepair
from mappa.domain import AAMProblem, MolecularEndpoint
from mappa.native_search import find_islands_native
from mappa.search_graph import frozen_value
from mappa.shared_seed_search import SharedSeedSearch


def cut_sessions(problem, config):
    condition = AdaptiveSearchCondition.uncut(problem, config)
    repair = FragmentRepair(condition.source, condition.target, condition.orbits)
    for cut in cut_sweep_items(problem.reactant.wbo, config.cut_floor):
        view = repair.for_cut(tuple(e for e in cut if condition.source.has_edge(*e)))
        yield SharedSeedSearch(problem, config, condition=replace(condition,
            source=view.source, cuts=tuple(cut), growth_replay=view))


def state_key(state):
    return state.mapping, state.islands, state.deferred_edges


def contents(graph):
    # Ignore duplicated histories, node IDs and seed-provenance multiplicity.
    # Retain exact correlated fragment records, constraints and stop reasons.
    states = {state_key(s) for s in graph.states}
    edges = {(state_key(graph.states[e.source]), state_key(graph.states[e.target]),
              frozen_value(e.match), e.preserved_bonds) for e in graph.transitions if e.match is not None}
    stops = {(state_key(graph.states[s.state]), s.reason, s.stage) for s in graph.stops}
    return states, edges, stops


@pytest.mark.parametrize('tol', [.5, 1.])
@pytest.mark.parametrize('cap', [1, 3, 100])
def test_exact_union_of_original_policies_including_frontier_caps(tol, cap):
    rng = np.random.default_rng(82197)
    for trial in range(6):
        endpoints = []
        for side in range(2):
            n = 7 + (side if trial % 2 else 0)
            elements = tuple(['C']*4 + ['O']*2 + ['H']*(n-6))
            w = np.triu(rng.choice([0., 0., 1., 1.5, 2.], size=(n,n)), 1)
            endpoints.append(MolecularEndpoint(elements, np.zeros((n,3)), w+w.T, str(side)))
        problem = AAMProblem(*endpoints)
        config = AAMSearchConfig(seed_count=3, branch_limit=cap, iso_tolerance=tol)
        sessions = list(cut_sessions(problem, config))
        child = sessions[trial % len(sessions)]
        while child.advance():
            pass
        actual = contents(child.builder.finish())
        reference = [set(), set(), set()]
        for order in child.orders:
            graph = find_islands_native(child.source, child.target, order,
                graph_floor=config.graph_floor, iso_tol=tol, max_branches=cap,
                p_orbits=child.condition.orbits, r_orbits=child.source_orbits,
                cuts=child.condition.cuts)
            for target, values in zip(reference, contents(graph)):
                target.update(values)
        assert actual == tuple(reference)
        assert nx.is_directed_acyclic_graph(nx.DiGraph((e.source,e.target)
            for e in child.builder.transitions))
        assert child.growth_calls == len(child.decisions)
        assert child.work == child.growth_calls + child.reused_states


def test_shared_growth_cache_and_saved_snapshot_preserve_groups():
    from test_fragment_choices import problem
    sessions = list(cut_sessions(problem(6), AAMSearchConfig(seed_count=3, iso_tolerance=1.)))
    for child in sessions:
        for _ in range(25):
            child.advance()
        prefix = child.snapshot()
        record = prefix.aam.graph.to_record()
        while child.advance():
            pass
        assert prefix.aam.graph.to_record() == record
        plain = SharedSeedSearch(child.problem, child.config,
            condition=replace(child.condition, growth_replay=None))
        while plain.advance():
            pass
        assert plain.snapshot().aam.graph == child.snapshot().aam.graph
    assert sum(c.reused_states for c in sessions) > 0


@pytest.mark.parametrize('cap', [1, 3, 100])
def test_element_exhaustion_skips_only_impossible_growth(cap, monkeypatch):
    """Other elements can remain available when a particular seed is impossible."""
    source = MolecularEndpoint(('C',)*5+('O',)*3, np.zeros((8,3)), np.zeros((8,8)))
    target = MolecularEndpoint(('C','O','N'), np.zeros((3,3)), np.zeros((3,3)))
    problem = AAMProblem(source,target)
    config = AAMSearchConfig(seed_count=3, branch_limit=cap, iso_tolerance=1.)
    pruned = next(cut_sessions(problem, config))
    while pruned.advance():
        pass
    reference = SharedSeedSearch(problem,config,condition=pruned.condition)
    monkeypatch.setattr(reference,'_can_seed',lambda branch,seed: seed not in branch.mapping)
    while reference.advance():
        pass
    assert contents(pruned.builder.finish()) == contents(reference.builder.finish())
    assert pruned.growth_calls < reference.growth_calls


@pytest.mark.parametrize('workers', [1, 2])
def test_cut_worker_backend_preserves_full_fragment_content_and_checkpoint_identity(tmp_path, workers):
    from mappa import search_aam
    from test_fragment_choices import problem
    config = AAMSearchConfig(seed_count=3, iso_tolerance=1., branch_limit=100)
    original = search_aam(problem(6), config, workers=workers, execution='reused_native')
    directory = tmp_path/'cuts'
    shared = search_aam(problem(6), config, workers=workers, execution='shared_policies',
                        intermediate_dir=directory, archive_format='checkpoint')
    # Conditioned symmetry is finalized too: compare complete records, not only scores.
    assert contents(shared.graph) == contents(original.graph)
    assert len(list(directory.glob('*.raw.pkl.gz'))) == shared.metrics.cut_count
    resumed = search_aam(problem(6), config, workers=workers, execution='shared_policies',
                         intermediate_dir=directory, archive_format='checkpoint', resume=True)
    assert resumed.graph == shared.graph
    with pytest.raises(ValueError, match='configuration differs'):
        search_aam(problem(6), config, workers=workers, execution='reused_native',
                   intermediate_dir=directory, archive_format='checkpoint', resume=True)
