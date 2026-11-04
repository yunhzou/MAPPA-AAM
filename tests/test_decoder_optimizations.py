"""Exactness checks for performance reductions, without permutation enumeration."""
from itertools import combinations
import time
from types import SimpleNamespace

import numpy as np
import pytest
import z3

from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.core.symmetry.event_counts import _event_model, _event_count_constraint
from mappa.mappa_core.core.symmetry.family import compile_path
from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
from mappa.mappa_core.post_processing.final_branches import FinalFamily
from mappa.mappa_core.post_processing.decoder.events import classify_bonds
from mappa.mappa_core.utils.bonds import is_metal_element
from reference_event_model import _event_model as reference_event_model
from test_event_patterns import problem
from growth_cases import saved_fragment


def test_composed_local_relations_preserve_injectivity():
    value = problem('OOOOO', [], [])
    first = saved_fragment(value, source_atoms=(0,1,2), target_atoms=(0,1,2))
    last = saved_fragment(value, source_atoms=(2,3,4), target_atoms=(2,3,4), anchors=((2,2),))
    family = FinalFamily(tuple(enumerate(range(5))), (),
                         (('relation',last), ('relation',first)), (.5,1.))
    compiled = compile_path(family.as_path(value), value, {}, source_atoms=(), complete_reference=False)
    images = [v[0] for v in compiled.values]
    # No global Distinct premise: injectivity follows from local joint choices.
    assert not any(z3.is_distinct(a) and a.num_args() == 5 for a in compiled.solver.assertions())
    compiled.solver.add(z3.Not(z3.Distinct(*images)))
    compiled.solver.set(timeout=2000)
    assert compiled.solver.check() == z3.unsat


def test_cached_pair_tables_and_boolean_counts_equal_scalar_tables():
    rng = np.random.default_rng(24)
    weights = [0., .3, .5, np.nextafter(.5, 0), np.nextafter(.5, 1), 1., 1.4]
    for elements in ('OOOO', 'VOOO'):
        edges = list(combinations(range(4), 2))
        value = problem(elements, [(a,b,float(rng.choice(weights))) for a,b in edges],
                        [(a,b,float(rng.choice(weights))) for a,b in edges])
        index = SignedEventIndex(value)
        constraints = FiniteConstraints()
        physical = [constraints.symbol(j for j, element in enumerate(elements)
                                       if element == elements[i]) for i in range(4)]
        constraints.distinct(tuple(value for value, _ in physical))
        solver = z3.Solver()
        constraints.emit(solver)
        values = [(constraints.expression(value), domain) for value, domain in physical]
        compiled = SimpleNamespace(values=values, physical_values=physical,
            constraints=constraints, solver=solver, representative=dict(enumerate(range(4))))
        old, total, metrics = reference_event_model(compiled, index, time.perf_counter()+5)
        for _ in range(2):  # Cold build followed by cache reuse.
            new, new_total, new_metrics = _event_model(compiled, index, time.perf_counter()+5)
            assert new_metrics['variable_pairs'] == metrics['variable_pairs']
            assert new_metrics['event_lower_bound'] == metrics['event_lower_bound']
            assert new_metrics['pair_table_entries'] <= metrics['pair_table_entries']
            differences = [total != new_total] + [old[k] != new.get(k, 0) for k in old]
            for count in range(8):
                differences.extend([
                    (total == count) != _event_count_constraint(new, count),
                    (total <= count) != _event_count_constraint(new, count, at_most=True)])
            compiled.solver.push(); compiled.solver.add(z3.Or(*differences))
            compiled.solver.set(timeout=2000)
            assert compiled.solver.check() == z3.unsat
            compiled.solver.pop()
        assert index._pair_tables


def test_sparse_scalar_validation_matches_dense_mapped_boundary_rule():
    rng = np.random.default_rng(718)
    elements = tuple('COHVOO')
    weights = [0., .3, np.nextafter(.3, 0), .5, np.nextafter(.5, 1), 1.]
    edges = list(combinations(range(6), 2))
    for _ in range(12):
        matrices = []
        for side in range(2):
            w = np.zeros((6,6))
            for a,b in edges:w[a,b]=w[b,a]=rng.choice(weights)
            matrices.append(w)
        r,p = matrices
        # Fixed bijection with partial variants; no mapping search/enumeration.
        mapping = {a:b for a,b in enumerate((4,1,5,3,0,2)) if rng.random() > .2}
        inv = {b:a for a,b in mapping.items()}
        for metal in (None, .3):
            def events(source, target, m):
                out=[]
                for a,b in edges:
                    if a not in m and b not in m:continue
                    tau = metal if metal is not None and any(is_metal_element(elements[i]) for i in (a,b)) else .5
                    w=source[a,b]
                    if w <= tau:continue
                    mapped = a in m and b in m
                    other=target[m[a],m[b]] if mapped else None
                    if not mapped or w-other > tau:out.append((a,b,float(w),None if other is None else float(other)))
                return out
            broken=events(r,p,mapping)
            formed=[(a,b,other,w) for a,b,w,other in events(p,r,inv)]
            expected=(broken,formed,sorted({a for row in broken for a in row[:2]}),
                      sorted({a for row in formed for a in row[:2]}))
            assert classify_bonds(mapping,r,p,elements_R=elements,elements_P=elements,
                                  metal_dwbo_threshold=metal)==expected


def test_empty_endpoints_keep_empty_scalar_event_result():
    w = np.zeros((0,0))
    assert classify_bonds({},w,w,elements_R=(),elements_P=(),metal_dwbo_threshold=.3) == ([],[],[],[])


def test_fragment_contributions_match_dense_events_for_all_small_partial_maps():
    from itertools import permutations
    from mappa.mappa_core.post_processing.hierarchical_mapping import FragmentEventCounter
    from test_partial_event_decoding import endpoint
    from mappa import AAMProblem
    rng = np.random.default_rng(719)
    weights = [0., .3, np.nextafter(.3, 1), .5, np.nextafter(.5, 1), 1., 1.8]
    for elements, target_elements in [('CCCC', 'CCCCC'), ('CoCCC', 'CoCCCC')]:
        elements = ('Co', 'C', 'C', 'C') if elements.startswith('Co') else tuple(elements)
        target_elements = ('Co', 'C', 'C', 'C', 'C') if target_elements.startswith('Co') else tuple(target_elements)
        def molecule(es):
            return endpoint(es, [(a,b,float(rng.choice(weights))) for a,b in combinations(range(len(es)),2)])
        value = AAMProblem(molecule(elements), molecule(target_elements))
        for value in (value, AAMProblem(value.product, value.reactant)):
            index = SignedEventIndex(value)
            counter = FragmentEventCounter(index)
            for size in range(4):
                for atoms in combinations(range(value.source_atom_count), size):
                    for images in permutations(range(value.target_atom_count), size):
                        if any(value.reactant.elements[a] != value.product.elements[b] for a,b in zip(atoms,images)):
                            continue
                        mapping = tuple(zip(atoms, images))
                        # Independent dense endpoint oracle, including unmapped boundaries.
                        expected = []
                        for ep, other, assignment in ((value.reactant, value.product, dict(mapping)),
                                (value.product, value.reactant, {b:a for a,b in mapping})):
                            pairs = []
                            for a,b in combinations(range(ep.atom_count), 2):
                                if a not in assignment and b not in assignment:
                                    continue
                                tau = .3 if any(is_metal_element(ep.elements[i]) for i in (a,b)) else .5
                                opposite = other.wbo[assignment[a],assignment[b]] if a in assignment and b in assignment else 0.
                                if ep.wbo[a,b] - opposite > tau:
                                    pairs.append((a,b))
                            expected.append(tuple(pairs))
                        for owners in (tuple((a,0) for a in atoms), tuple((a,a) for a in atoms),
                                       tuple((a,a % 2) for a in atoms)):
                            total, internal = counter.count(mapping, owners)
                            assert counter.events(mapping, owners) == tuple(expected)
                            assert total == sum(map(len, expected))
                            r_owner = dict(owners)
                            p_owner = {p:r_owner[r] for r,p in mapping}
                            assert internal == sum(a in own and b in own and own[a] == own[b]
                                for events, own in zip(expected, (r_owner,p_owner)) for a,b in events)
                            before = counter.bond_comparisons
                            assert counter.count(mapping, owners) == (total, internal)
                            assert counter.bond_comparisons == before


def test_decoder_uses_cached_events_for_candidate_identity(monkeypatch):
    import mappa.mappa_core.post_processing.decoder.events as module
    from mappa.mappa_core.post_processing.hierarchical_mapping import refine_mapping
    from test_fragment_conversion import saved, endpoint
    ep = endpoint('COHH', [(0,1,1), (1,2,1), (1,3,1)])
    other = endpoint('COHH', [(0,2,1), (1,2,1), (1,3,1)])
    source = saved(ep, other, {i:i for i in range(4)})
    # The independent public description remains a separate correctness oracle.
    expected = SignedEventIndex(source.problem).describe(list(range(4)))
    def forbidden(*args, **kwargs):
        pytest.fail('decoder repeated whole-mapping bond counting')
    import pytest
    monkeypatch.setattr(module, 'classify_bonds', forbidden)
    monkeypatch.setattr(module.SignedEventIndex, 'counts', forbidden)
    result = refine_mapping(source)
    assert result.candidates[0].id == expected['id']
    assert result.candidates[0].endpoint_events == expected['endpoint_events']
    assert result.candidates[0].total == expected['total']
    assert result.counting_statistics['branch_counts'] == len(result.family_reports)


def test_many_combinations_reuse_fragment_and_boundary_counts():
    from itertools import product
    from mappa.mappa_core.post_processing.hierarchical_mapping import FragmentEventCounter
    import json
    # Twelve two-atom fragments, each with two internal correspondences.
    # Only neighboring fragments interact, so each link needs four evaluations.
    count = 12
    edges_r = [(2*i, 2*i+1, 1.) for i in range(count)]
    edges_r += [(2*i+1, 2*i+2, 1.) for i in range(count-1)]
    edges_p = [(2*i, 2*i+1, 1.) for i in range(count)]
    edges_p += [(2*i, 2*i+2, 1.) for i in range(count-1)]
    value = problem('C'*(2*count), edges_r, edges_p)
    owners = tuple((a, a//2) for a in range(2*count))
    mappings = [tuple((2*i+a, 2*i+(a ^ choice)) for i,choice in enumerate(choices) for a in (0,1))
                for choices in product((0,1), repeat=count)]
    start = time.process_time()
    direct = []
    for mapping in mappings:
        broken, formed, _, _ = classify_bonds(dict(mapping), value.reactant.wbo, value.product.wbo,
            elements_R=value.reactant.elements, elements_P=value.product.elements, metal_dwbo_threshold=.3)
        direct.append(len(broken)+len(formed))
    direct_cpu = time.process_time() - start
    start = time.process_time()
    counter = FragmentEventCounter(SignedEventIndex(value))
    actual = [counter.count(mapping, owners)[0] for mapping in mappings]
    cached_cpu = time.process_time() - start
    assert actual == direct
    assert len(counter.boundaries) == 4*(count-1)
    assert counter.bond_comparisons == 2*(2*count) + 2*4*(count-1)
    assert len(counter.topologies) == 1
    print(json.dumps(dict(combinations=len(mappings), direct_cpu=direct_cpu,
        cached_cpu=cached_cpu, statistics=counter.statistics)), flush=True)


@pytest.mark.parametrize('edges,tolerance,expected_counts', [
    ([(0, 1, 1.), (1, 2, 1.), (2, 3, 1.8), (3, 0, 1.8)], 1., (0, 2, 4)),
    ([(0, 1, 1.), (0, 2, 1.6), (0, 3, 2.2)], .7, (0, 2, 3)),
])
def test_local_spectrum_preserves_joint_count_gaps_and_cut_endpoint_correlation(
        edges, tolerance, expected_counts):
    from itertools import product
    from mappa.mappa_core.core.symmetry.event_counts import LocalEventCounts
    from test_relaxed_isomorphism import self_oracle

    value = problem('CCCC', edges, edges)
    fragment = saved_fragment(value, domains=tuple((a, (a,)) for a in range(4)),
        accepted_edges=tuple((a, b) for a, b, _ in edges), self_tolerance=tolerance)
    expected = {}
    for source, target in product(self_oracle(fragment.source_self_automorph),
                                  self_oracle(fragment.target_self_automorph)):
        images = tuple(target[source.index(atom)] for atom in range(4))
        total = sum(abs(value.reactant.wbo[a, b] - value.product.wbo[images[a], images[b]]) > .5
                    for a, b in combinations(range(4), 2))
        expected.setdefault((images[0],), set()).add(int(total))
    counter = LocalEventCounts(SignedEventIndex(value))
    actions = (('relation', fragment),)
    result = counter.spectrum(fragment.witness, actions, endpoints=(0,))
    assert result.minimum == 0
    assert result.counts == expected_counts
    assert dict(result.endpoint_counts) == {key: tuple(sorted(counts)) for key, counts in expected.items()}
    before = counter.solver_queries
    assert counter.spectrum(fragment.witness, actions, endpoints=(0,)) is result
    assert counter.solver_queries == before and counter.cache_hits == 1
    assert counter.spectrum(fragment.witness, actions, allowance=0).counts == (0,)
    assert counter.solver_queries == before
    limited = counter.spectrum(fragment.witness, actions, endpoints=(0,), allowance=2)
    assert limited.counts == (0, 2)
    assert dict(limited.endpoint_counts) == {
        key: tuple(sorted(counts & {0, 1, 2})) for key, counts in expected.items()
        if counts & {0, 1, 2}}


def test_local_spectrum_excludes_boundary_events_and_uses_original_placement_weights():
    from mappa.mappa_core.core.symmetry.event_counts import LocalEventCounts

    matching = problem('CCCCO', [(0, 1, 1.), (2, 3, 1.), (0, 4, 1.)],
                       [(0, 1, 1.), (2, 3, 1.)])
    fragment = saved_fragment(matching, source_atoms=(0, 1), accepted_edges=((0, 1),))
    original = problem('CCCCO', [(0, 1, 1.), (2, 3, 1.), (0, 4, 1.)],
                       [(0, 1, 1.), (2, 3, 1.8)])
    actions = (('relation', fragment), ('group', ((2, 3, 0, 1, 4),)))
    counter = LocalEventCounts(SignedEventIndex(original))
    result = counter.spectrum(fragment.witness, actions, endpoints=(0,))
    assert result.minimum == 0 and result.counts == (0, 1)
    assert dict(result.endpoint_counts) == {(0,): (0,), (1,): (0,), (2,): (1,), (3,): (1,)}
    # The source boundary to O is deliberately outside this internal response.
    assert counter.spectrum(fragment.witness, actions[:1]).counts == (0,)


def test_local_spectrum_keeps_overlapping_actions_ordered_and_endpoints_joint():
    from itertools import product
    from mappa.mappa_core.core.symmetry.event_counts import LocalEventCounts

    value = problem('CCC', [(0, 1, 1.), (1, 2, 1.8)], [(0, 1, 1.), (1, 2, 1.8)])
    fragments = tuple(saved_fragment(value, source_atoms=atoms,
        domains=tuple((atom, (atom,)) for atom in atoms), accepted_edges=(atoms,))
        for atoms in ((0, 1), (1, 2)))
    counter = LocalEventCounts(SignedEventIndex(value))
    for order in (fragments, fragments[::-1]):
        expected = {}
        for swaps in product((False, True), repeat=2):
            images = (0, 1)
            for fragment, exchange in zip(order, swaps):
                if exchange:
                    left, right = fragment.source_atoms
                    images = tuple(right if image == left else left if image == right else image
                                   for image in images)
            count = int(abs(value.reactant.wbo[0, 1] - value.product.wbo[images]) > .5)
            expected[images] = (count,)
        result = counter.spectrum(((0, 0), (1, 1)),
            tuple(('relation', fragment) for fragment in order), endpoints=(0, 1))
        assert dict(result.endpoint_counts) == expected
        assert all(left != right for left, right in dict(result.endpoint_counts))
    assert len(counter.minima) == 2


@pytest.mark.parametrize('endpoints', [(), (0,)])
def test_local_spectrum_reuses_widest_count_and_endpoint_support(endpoints):
    from mappa.mappa_core.core.symmetry.event_counts import LocalEventCounts

    edges = [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.8), (3, 0, 1.8)]
    value = problem('CCCC', edges, edges)
    fragment = saved_fragment(value, domains=tuple((a, (a,)) for a in range(4)),
        accepted_edges=tuple((a, b) for a, b, _ in edges))
    actions = (('relation', fragment),)
    expanded = LocalEventCounts(SignedEventIndex(value))
    narrow = expanded.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=2)
    assert narrow.counts == (0, 2)
    before = expanded.solver_queries
    full = expanded.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=None)
    assert full.counts == (0, 2, 4)
    assert expanded.solver_queries > before
    before = expanded.solver_queries
    assert expanded.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=2) == narrow
    minimum = expanded.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=0)
    assert minimum.counts == (0,)
    assert all(counts == (0,) for _, counts in minimum.endpoint_counts)
    assert expanded.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=4) is full
    assert expanded.solver_queries == before
    assert len(expanded.cache) == 1

    wide_first = LocalEventCounts(SignedEventIndex(value))
    wide_first.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=4)
    before = wide_first.solver_queries
    assert wide_first.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=2) == narrow
    assert wide_first.spectrum(fragment.witness, actions, endpoints=endpoints, allowance=0) == minimum
    assert wide_first.solver_queries == before
    assert len(wide_first.cache) == 1


@pytest.mark.parametrize('element,threshold', [('O', .5), ('Fe', .3)])
@pytest.mark.parametrize('above', [False, True])
def test_local_spectrum_respects_strict_ordinary_and_metal_thresholds(element, threshold, above):
    from mappa.mappa_core.core.symmetry.event_counts import LocalEventCounts

    weight = np.nextafter(threshold, np.inf) if above else threshold
    value = problem(('C', element), [(0, 1, weight)], [])
    result = LocalEventCounts(SignedEventIndex(value)).spectrum(((0, 0), (1, 1)), ())
    assert result.minimum == int(above)
    assert result.counts == (int(above),)


def test_local_spectrum_reports_unknown_and_rejects_inconsistent_relations(monkeypatch):
    from mappa.mappa_core.core.symmetry.event_counts import LocalEventCounts

    value = problem('CC', [(0, 1, 1.)], [])
    counter = LocalEventCounts(SignedEventIndex(value))
    with pytest.raises(ValueError, match='inconsistent'):
        counter.spectrum(((0, 0), (1, 0)), ())
    monkeypatch.setattr(z3.Solver, 'check', lambda self: z3.unknown)
    monkeypatch.setattr(z3.Solver, 'reason_unknown', lambda self: 'test exhaustion')
    with pytest.raises(RuntimeError, match='inconclusive.*test exhaustion'):
        counter.spectrum(((0, 0), (1, 1)), ())


@pytest.mark.parametrize('source_weight,expected_placed', [(1., 0), (.4, 1)])
def test_sweep_placement_only_cancels_events_present_in_fixed_boundary_ledger(
        source_weight, expected_placed):
    from mappa import AAMProblem, decode_events
    from test_event_interpretation import endpoint, saved

    source = endpoint('CC', ((0, 1, 1.),))
    target = endpoint('CCCC', ((0, 1, 1.), (2, 3, 1.)))
    matching = AAMProblem(source, target)
    fragments = tuple(saved_fragment(matching, source_atoms=(atom,), target_atoms=(atom,),
                                    witness={atom: atom}) for atom in (0, 1))
    catalogue = saved(source, target, {0: 0, 1: 1}, fragments=fragments,
                      pairs=(((0,), (0,)), ((1,), (1,))))
    catalogue.families[0].actions += (('group', ((2, 3, 0, 1),)),)
    catalogue.families[0].provenance = [{'cuts': ((0, 1),)}]
    original = AAMProblem(endpoint('CC', ((0, 1, source_weight),)),
                          endpoint('CCCC', ((0, 1, .4), (2, 3, 1.))))
    decoded = decode_events(catalogue, event_problem=original)
    assert decoded.minimum_events == 0
    assert decoded.attainable_counts(0) == (0, 1)
    assert decoded.tables(0).boundary_count == int(source_weight > .5)
    choices = (None, None, (1,))
    assert decoded.realize(0, choices) == {0: 2, 1: 3}
    report = decoded.contributions(0, choices, include_events=True)
    assert report['total'] == expected_placed == len(report['events'])
    assert min(report['counts'].values()) >= 0
    assert sum(next(decoded.counts(0, (choices,)))) == expected_placed


def test_independent_bins_do_not_claim_a_joint_mapping_witness_for_global_placements():
    from mappa import AAMProblem, decode_events
    from test_event_interpretation import endpoint, saved

    elements = ('C', 'O', 'O', 'O') * 2
    bonds = tuple((base, base + i, 1.) for base in (0, 4) for i in (1, 2, 3))
    source = endpoint(elements, bonds + ((0, 4, 1.),))
    target = endpoint(elements, bonds)
    matching = AAMProblem(source, target)
    fragments = tuple(saved_fragment(matching, source_atoms=tuple(range(base, base + 4)),
        target_atoms=tuple(range(base, base + 4)),
        witness={atom: atom for atom in range(base, base + 4)},
        accepted_edges=tuple((base, base + i) for i in (1, 2, 3))) for base in (0, 4))
    catalogue = saved(source, target, dict(enumerate(range(8))), fragments=fragments,
        pairs=(((0, 1, 2, 3), (0, 1, 2, 3)), ((4, 5, 6, 7), (4, 5, 6, 7))))
    catalogue.families[0].actions += (('group', ((4, 5, 6, 7, 0, 1, 2, 3),)),)
    catalogue.families[0].provenance = [{'cuts': ((0, 4),)}]
    original = AAMProblem(source, endpoint(elements,
        tuple((a, b, 1. if a < 4 else 1.8) for a, b, _ in bonds)))
    for slack in (None, 6):
        decoded = decode_events(catalogue, event_problem=original, event_slack=slack)
        assert decoded.minimum_events == 1
        assert decoded.attainable_counts(0) == (1, 4, 7)
        for word in ((), (1,)):
            choice = (None, None, word)
            mapping = decoded.realize(0, choice)
            assert len(mapping) == len(set(mapping.values())) == 8
            # One of the two target stars always has three changed bonds.
            assert decoded.contributions(0, choice)['total'] == 4
            assert sorted(mapping.values()) == list(range(8))


def test_sweep_cancellation_uses_original_target_boundary_owner_region():
    from mappa import AAMProblem, decode_events
    from test_event_interpretation import endpoint, saved

    edges = ((0, 1, 1.), (0, 2, 1.), (1, 2, 1.))
    ep = endpoint('CCC', edges)
    matching = AAMProblem(ep, ep)
    fragments = tuple(saved_fragment(matching, source_atoms=(atom,), target_atoms=(atom,),
                                    witness={atom: atom}) for atom in range(3))
    catalogue = saved(ep, ep, dict(enumerate(range(3))), fragments=fragments,
        pairs=tuple(((atom,), (atom,)) for atom in range(3)))
    catalogue.families[0].actions += (('group', ((1, 2, 0),)),)
    catalogue.families[0].provenance = [{'cuts': ((0, 1),)}]
    original = AAMProblem(endpoint('CCC', ((0, 1, 1.), (0, 2, .4), (1, 2, .4))),
                          endpoint('CCC', ((0, 1, .4), (0, 2, .4), (1, 2, 1.))))
    decoded = decode_events(catalogue, event_problem=original)
    assert decoded.attainable_counts(0) == (0, 2)
    report = decoded.contributions(0, (None, None, None, (1,)), include_events=True)
    assert report['total'] == 0
    assert report['boundary'] == {}
    assert report['events'] == ()
