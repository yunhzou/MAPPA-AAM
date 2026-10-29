"""Exact target actions must transport future atoms, not discard their matches."""
from itertools import product, permutations, combinations

import numpy as np
import pytest

from mappa import FragmentMatchConfig, match_fragment
from mappa.mappa_core.core.matching.state import _SymBlock, _SymCand
from mappa.mappa_core.core.matching.dedupe import _dedup_sym_cands
from mappa.mappa_core.core.matching.extend import _extend_sym_cands
from mappa.mappa_core.core.matching.orbits import _nauty_orbits
from mappa.mappa_core.core.matching.state import candidate_from_record
from mappa.mappa_core.core.matching.symmetry import SymmetryWorkspace
from test_growth_compression_exactness import graph, represented, assert_prefixes


def rings():
    cycle = [(i, (i + 1) % 6, 1) for i in range(6)]
    source = graph(['C'] * 6 + ['H'] * 6, cycle + [(i, 6 + i, 1) for i in range(6)])
    target = graph(['C'] * 6 + ['H'] * 12,
                   cycle + [(i, 6 + 2 * i + h, 1) for i in range(6) for h in range(2)])
    return source, target


def ring_mappings():
    # Independent construction: six rotations, two directions, two H per C.
    result = set()
    for offset, direction in product(range(6), (-1, 1)):
        carbons = [(offset + direction * i) % 6 for i in range(6)]
        for choices in product(range(2), repeat=6):
            images = carbons + [6 + 2 * p + h for p, h in zip(carbons, choices)]
            result.add(tuple(enumerate(images)))
    return result


@pytest.mark.parametrize('order', [list(range(12)), [a for i in range(6) for a in (i, 6 + i)]])
def test_ring_every_prefix_preserves_all_768_matches(order):
    source, target = rings()
    ro, po = _nauty_orbits(source, wbo_tol=1), _nauty_orbits(target, wbo_tol=1)
    candidates = [_SymCand({0: 0}, (_SymBlock((0,), tuple(range(6)), False),))]
    fragment = {0}
    for atom in order[1:]:
        candidates = _extend_sym_cands(candidates, fragment, atom, source, target,
                                       {}, 1, None, r_orbits=ro, p_orbits=po)
        fragment.add(atom)
        expected = {tuple((r, p) for r, p in mapping if r in fragment) for mapping in ring_mappings()}
        assert represented(candidates, target, 1) == expected
    assert len(candidates) == 1


@pytest.mark.parametrize('native', ['0', '1'])
@pytest.mark.smoke
def test_public_ring_match_is_one_branch_with_768_assignments(native, monkeypatch):
    monkeypatch.setenv('MAPPA_NATIVE', native)
    source, target = rings()
    result = match_fragment(source, target, seed=0,
                            config=FragmentMatchConfig(iso_tolerance=1, branch_limit=100))
    assert not result.capped
    assert len(result.matches) == 1
    candidate = candidate_from_record(result.matches[0].symmetry)
    assert set(candidate.mapping) == set(source)
    assert represented([candidate], target, 1) == ring_mappings()


def test_target_actions_fix_locked_images():
    _, target = rings()
    candidate = _SymCand({0: 0})
    generators = SymmetryWorkspace(target, 1).conditioned_generators(candidate, ((99, 2), (100, 10)))
    assert all(g[2] == 2 and g[10] == 10 for g in generators)


def test_weak_weight_breaks_apparent_ring_symmetry():
    source, target = rings()
    target.graph['wbo_matrix'][0, 2] = target.graph['wbo_matrix'][2, 0] = .05
    generators = SymmetryWorkspace(target, 1).conditioned_generators(_SymCand({0: 0}), ())
    weights = target.graph['wbo_matrix']
    assert all(np.array_equal(weights, weights[np.ix_(g, g)]) for g in generators)


def test_new_branch_is_not_constrained_by_old_automorph_domain_metadata():
    target = graph(['C'] * 4, [(i, (i + 1) % 4, 1) for i in range(4)])
    first = _SymCand({0: 0, 1: 1}, automorph_blocks=(_SymBlock((0,), (0, 2), False),))
    second = _SymCand({0: 1, 1: 2})
    candidates = _dedup_sym_cands([first, second], target, target)
    assert len(candidates) == 1


def test_asymmetric_future_neighbors_keep_non_equivalent_matches():
    source = graph(['C', 'C', 'O'], [(0, 1, 1), (1, 2, 1)])
    target = graph(['C', 'C', 'C', 'O'], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    assert_prefixes(source, target, [0, 1, 2])


def test_decoding_transports_later_fragments_with_the_earlier_group():
    from mappa import AAMProblem, query_path
    from mappa.alignment.branch import find_islands
    from mappa.search_symmetry import finalize_graph_symmetry
    from test_family_query import endpoint

    source = graph(['C'] * 4, [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    target = graph(['C'] * 4, [(i, (i + 1) % 4, 1) for i in range(4)])
    problem = AAMProblem(endpoint('CCCC', source.edges), endpoint('CCCC', target.edges))
    found = find_islands(source, target, [0, 1, 2, 3], iso_tol=.2, max_branches=100)
    found, _ = finalize_graph_symmetry(found, target, iso_tolerance=.2)
    paths = list(found.paths())
    assert paths and all(len(path.fragments) > 1 for path in paths)
    for path in paths:
        first = path.hierarchy.fragments[0]
        outside = set(target) - {p for _, p in first.representative_assignments}
        assert any(g.images[p] != p for g in first.target_generators for p in outside)
        required = set()
        for fragment in path.hierarchy.fragments:
            deferred = {tuple(sorted(e)) for e in fragment.deferred_edges}
            required.update((a, b) for a, b in combinations(sorted(fragment.r_atoms), 2)
                            if source.has_edge(a, b) and (a, b) not in deferred)
        for images in permutations(range(4)):
            expected = all(target.has_edge(images[a], images[b]) for a, b in required)
            status, _ = query_path(path, problem, dict(enumerate(images)),
                                   source_atoms=tuple(range(4)), timeout_ms=2000)
            assert status == ('recovered' if expected else 'not_recovered')


def test_locked_ring_then_free_hydrogen():
    """Run the real scheduler and enumerate this small saved decoder relation."""
    import json
    import time
    import z3
    from mappa import AAMProblem
    from mappa.alignment.branch import find_islands
    from mappa.family_query import compile_path
    from mappa.search_symmetry import finalize_graph_symmetry
    from test_family_query import endpoint

    ring, target = rings()
    source = graph(['C'] * 6 + ['H'] * 7, [(a, b, 1) for a, b in ring.edges])
    problem = AAMProblem(endpoint('C' * 6 + 'H' * 7, source.edges),
                         endpoint('C' * 6 + 'H' * 12, target.edges))
    start = time.process_time()
    found = find_islands(source, target, list(range(13)), iso_tol=1, max_branches=100)
    found, _ = finalize_graph_symmetry(found, target, iso_tolerance=1)
    paths = list(found.paths())
    assert len(paths) == 1
    path = paths[0]
    assert [set(f['fragment']) for f in path.fragments] == [set(range(12)), {12}]
    first, second = path.hierarchy.fragments
    locked = dict(first.representative_assignments)
    remaining = set(range(6, 18)) - set(locked.values())
    assert len(remaining) == 6
    assert len(second.symmetry_domains) == 1
    assert set(second.symmetry_domains[0].p_atoms) == remaining
    search_cpu = time.process_time() - start

    expected = {mapping + ((12, h),) for mapping in ring_mappings()
                for h in set(range(6, 18)) - {p for _, p in mapping}}
    assert len(expected) == 4608
    compiled = compile_path(path, problem, {}, source_atoms=tuple(range(13)),
                            complete_reference=False)
    compiled.solver.set(timeout=3000)
    decoded = set()
    decoding_start = time.process_time()
    while True:
        status = compiled.solver.check()
        assert status != z3.unknown, compiled.solver.reason_unknown()
        if status == z3.unsat:
            break
        model = compiled.solver.model()
        # Use the production realization code, including its independent
        # injectivity, element, and retained-bond validation.
        mapping = tuple(map(tuple, compiled.realize(model)['mapping']))
        assert mapping in expected and mapping not in decoded
        decoded.add(mapping)
        compiled.solver.add(z3.Or(*(value[0] != p
            for value, (_r, p) in zip(compiled.values, mapping))))
        if len(decoded) % 768 == 0:
            print(f'decoded {len(decoded)}', flush=True)
    assert decoded == expected
    result = dict(compressed_paths=len(paths), fragment_sizes=[12, 1],
                  locked_witness=locked, remaining_H=sorted(remaining),
                  H_pool=sorted(second.symmetry_domains[0].p_atoms),
                  decoded_mappings=len(decoded), extra_H_targets=sorted({m[-1][1] for m in decoded}),
                  choices_per_ring_mapping=6, search_cpu_seconds=search_cpu,
                  decoder_cpu_seconds=time.process_time()-decoding_start,
                  terminal_solver_status='unsat', exact_oracle_equality=True)
    print(json.dumps(result, sort_keys=True), flush=True)
