"""Independent oracles for fragment counts, fixed boundaries and one sweep cut."""
from dataclasses import replace
from itertools import combinations, permutations, product
import json

import numpy as np
import pytest

from mappa import AAMProblem, decode_events
from mappa.mappa_core.post_processing.decoder import EventDecodeConfig
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.utils.io.records import decoded_events_record, problem_from_record
from growth_cases import saved_fragment
from test_event_interpretation import endpoint, saved


def approximate_counts(problem, representative, partition, mapping, cut=()):
    """Scalar scientific rule, independent of response tables and native code."""
    source, target = problem.reactant, problem.product
    owners = {atom: owner for owner, atoms in enumerate(partition) for atom in atoms}
    target_owners = {representative[a]: owner for a, owner in owners.items()}
    counts = dict(breaking=0, forming=0, strengthening=0, weakening=0)

    def threshold(elements, a, b):
        return .3 if 'Fe' in (elements[a], elements[b]) or 'V' in (elements[a], elements[b]) else .5

    def internal(old, new, tau):
        if old - new > tau:
            counts['weakening' if new > tau else 'breaking'] += 1
        elif new - old > tau:
            counts['strengthening' if old > tau else 'forming'] += 1

    for atoms in partition:
        for a, b in combinations(atoms, 2):
            internal(source.wbo[a, b], target.wbo[mapping[a], mapping[b]],
                     threshold(source.elements, a, b))
    for ep, regions, category in ((source, owners, 'breaking'),
                                   (target, target_owners, 'forming')):
        for a, b in combinations(range(ep.atom_count), 2):
            left, right = regions.get(a, -1), regions.get(b, -1)
            if left != right and ep.wbo[a, b] > threshold(ep.elements, a, b):
                counts[category] += 1
    assert len(cut) <= 1
    for a, b in cut:
        if a not in mapping or b not in mapping or owners[a] == owners[b]:
            continue
        tau = threshold(source.elements, a, b)
        old, new = source.wbo[a, b], target.wbo[mapping[a], mapping[b]]
        counts['breaking'] -= int(old > tau)
        x, y = mapping[a], mapping[b]
        counted_target = (target_owners.get(x, -1) != target_owners.get(y, -1)
                          and target.wbo[x, y] > threshold(target.elements, x, y))
        counts['forming'] -= int(counted_target)
        internal(old, new, tau)
    return counts


def exact_count(problem, mapping):
    return sum(abs(problem.reactant.wbo[a, b] -
                   problem.product.wbo[mapping[a], mapping[b]]) > .5
               for a, b in combinations(mapping, 2))


def self_permutations(weights, tolerance=1., anchors=()):
    pairs = tuple(combinations(range(len(weights)), 2))
    return tuple(p for p in permutations(range(len(weights)))
        if all(p[a] == a for a in anchors)
        and all((weights[a, b] > .2) == (weights[p[a], p[b]] > .2)
                and abs(weights[a, b] - weights[p[a], p[b]]) <= tolerance for a, b in pairs))


def two_sided_choices(source, target=None, *, tolerance=1., anchors=()):
    target = source if target is None else target
    for alpha, beta in product(self_permutations(source, tolerance, anchors),
                                self_permutations(target, tolerance, anchors)):
        inverse = {image: atom for atom, image in enumerate(alpha)}
        mapping = {a: beta[inverse[a]] for a in range(len(source))}
        yield ({'source_permutation': alpha, 'target_permutation': beta},), mapping


@pytest.fixture
def ring_catalogue():
    ring = ((0, 1, 1.), (1, 2, 1.), (2, 3, 1.8), (3, 0, 1.8))
    ep = endpoint('CCCC', ring)
    problem = AAMProblem(ep, ep)
    edges = tuple(tuple(sorted((a, b))) for a, b, _ in ring)
    identity = dict(enumerate(range(4)))
    fragment = saved_fragment(problem, domains=tuple((a, (a,)) for a in range(4)),
                              accepted_edges=edges)
    catalogue = saved(ep, ep, identity, fragments=(fragment,))
    for offset in (2, 1):
        mapping = {a: (a + offset) % 4 for a in range(4)}
        fixed = saved_fragment(problem, witness=mapping,
                              anchors=tuple(mapping.items()), accepted_edges=edges)
        catalogue.families.append(saved(ep, ep, mapping, fragments=(fixed,)).families[0])
    next(iter(catalogue.branches.values())).update((1, 2))
    catalogue.path_count = 3
    return catalogue


def test_correlated_internal_bins_have_gaps_and_window_is_inclusive(ring_catalogue):
    catalogue = ring_catalogue
    choices = list(two_sided_choices(catalogue.problem.reactant.wbo))
    expected = [approximate_counts(catalogue.problem, dict(catalogue.families[0].mapping),
                 ((0, 1, 2, 3),), mapping) for _, mapping in choices]
    assert len(choices) == 64
    assert {sum(row.values()) for row in expected} == {0, 2, 4}
    unrestricted = decode_events(catalogue, event_slack=None)
    assert unrestricted.attainable_counts(0) == (0, 2, 4)
    decoded = decode_events(catalogue, event_slack=2)
    assert decoded.minimum_events == 0 and decoded.event_window == 2
    assert decoded.retained_family_ids == (0, 2)
    assert decoded.attainable_counts(0) == (0, 2)
    for (choice, mapping), counts in zip(choices, expected):
        if sum(counts.values()) <= 2:
            assert decoded.realize(0, choice) == mapping
            assert decoded.contributions(0, choice)['counts'] == counts
        else:
            with pytest.raises(ValueError, match='window'):
                decoded.realize(0, choice)
    assert decoded.contributions(2)['total'] == 2
    with pytest.raises(ValueError, match='window'):
        decoded.realize(1)
    default = decode_events(catalogue)
    assert EventDecodeConfig().event_slack == 6
    assert default.event_window == 6 and default.attainable_counts(0) == (0, 2, 4)
    batch = list(unrestricted.counts(0, (choice for choice, _ in choices), categories=True))
    assert [tuple(map(int, row)) for row in batch] == [tuple(row.values()) for row in expected]


def test_internal_excess_one_can_win_by_preserving_the_sweep_edge():
    source = endpoint('CCCO', ((0, 1, 1.), (0, 2, 1.4), (1, 2, 1.), (1, 3, 1.)))
    target = endpoint('CCCO', ((0, 1, 1.), (0, 2, 1.8), (1, 2, 1.), (2, 3, 1.)))
    problem = AAMProblem(source, target)
    fragment = saved_fragment(problem, source_atoms=(0, 1, 2), anchors=((0, 0),),
        domains=((0, (0,)), (1, (1,)), (2, (2,))), accepted_edges=((0, 1), (0, 2), (1, 2)))
    oxygen = saved_fragment(problem, source_atoms=(3,), anchors=((3, 3),))
    identity = dict(enumerate(range(4)))
    partition = ((0, 1, 2), (3,))
    catalogue = saved(source, target, identity, fragments=(fragment, oxygen),
        pairs=tuple((atoms, atoms) for atoms in partition))
    catalogue.families[0].provenance = [{'cuts': ((1, 3),)}]
    swapped = {0: 0, 1: 2, 2: 1, 3: 3}
    ordinary = approximate_counts(problem, identity, partition, swapped)
    corrected = approximate_counts(problem, identity, partition, swapped, ((1, 3),))
    assert sum(ordinary.values()) == 3 and sum(corrected.values()) == 1
    assert sum(approximate_counts(problem, identity, partition, identity, ((1, 3),)).values()) == 2
    unrestricted = decode_events(catalogue, event_slack=None)
    assert unrestricted.attainable_counts(0) == (1, 2)
    choice = ({'target_permutation': (0, 2, 1)}, None)
    batch = list(unrestricted.counts(0, (None, choice), categories=True))
    expected = [approximate_counts(problem, identity, partition, mapping, ((1, 3),))
                for mapping in (identity, swapped)]
    assert [tuple(map(int, row)) for row in batch] == [tuple(row.values()) for row in expected]
    decoded = decode_events(catalogue, event_slack=0)
    assert decoded.minimum_events == 1 and decoded.attainable_counts(0) == (1,)
    assert decoded.realize(0, choice) == swapped
    assert decoded.contributions(0, choice)['counts'] == corrected
    with pytest.raises(ValueError, match='window'):
        decoded.realize(0)


def test_ring_to_chain_counts_its_internal_cut_once():
    source = endpoint('CCCC', ((0, 1, 1.), (1, 2, 1.), (2, 3, 1.), (0, 3, 1.)))
    target = endpoint('CCCC', ((0, 1, 1.), (1, 2, 1.), (2, 3, 1.)))
    problem = AAMProblem(source, target)
    fragment = saved_fragment(problem, domains=tuple((a, (a,)) for a in range(4)),
                              accepted_edges=((0, 1), (1, 2), (2, 3)))
    catalogue = saved(source, target, dict(enumerate(range(4))), fragments=(fragment,))
    catalogue.families[0].provenance = [{'cuts': ((0, 3),)}]
    decoded = decode_events(catalogue, event_slack=0)
    assert decoded.minimum_events == 1 and decoded.attainable_counts(0) == (1,)
    for choice, mapping in two_sided_choices(source.wbo, target.wbo):
        assert exact_count(problem, mapping) == 1
        assert decoded.contributions(0, choice)['total'] == 1


def test_ordinary_ring_oxygen_boundary_uses_approximate_bin_for_reporting_too():
    ep = endpoint('CCCO', ((0, 1, 1.), (0, 2, 1.), (1, 2, 1.), (0, 3, 1.)))
    problem = AAMProblem(ep, ep)
    fragment = saved_fragment(problem, source_atoms=(0, 1, 2),
        domains=((0, (0,)), (1, (1,)), (2, (2,))), accepted_edges=((0, 1), (0, 2), (1, 2)))
    oxygen = saved_fragment(problem, source_atoms=(3,), anchors=((3, 3),))
    identity = dict(enumerate(range(4)))
    catalogue = saved(ep, ep, identity, fragments=(fragment, oxygen),
                      pairs=(((0, 1, 2), (0, 1, 2)), ((3,), (3,))))
    decoded = decode_events(catalogue)
    assert exact_count(problem, identity) == 0
    assert decoded.minimum_events == 2 and decoded.event_window == 8
    assert decoded.attainable_counts(0) == (2,)
    assert decoded.contributions(0)['counts'] == dict(breaking=1, forming=1, strengthening=0, weakening=0)
    assert decoded.evaluate(0)['total'] == 2
    assert sum(next(decoded.counts(0, [None]))) == 2


def test_distinct_cut_histories_never_accumulate_two_corrections():
    ep = endpoint('CCO', ((0, 1, 1.), (0, 2, 1.), (1, 2, 1.)))
    problem = AAMProblem(ep, ep)
    fragment = saved_fragment(problem, source_atoms=(0, 1), accepted_edges=((0, 1),))
    oxygen = saved_fragment(problem, source_atoms=(2,), anchors=((2, 2),))
    catalogue = saved(ep, ep, {0: 0, 1: 1, 2: 2}, fragments=(fragment, oxygen),
                      pairs=(((0, 1), (0, 1)), ((2,), (2,))))
    catalogue.families[0].provenance = [{'cuts': ((0, 2),)}, {'cuts': ((1, 2),)}]
    decoded = decode_events(catalogue, event_slack=0)
    assert set(decoded.inventory.cuts[0]) == {((0, 2),), ((1, 2),)}
    assert decoded.minimum_events == 2
    for partition in range(2):
        assert decoded.attainable_counts(0, partition=partition) == (2,)
        assert decoded.contributions(0, partition=partition)['total'] == 2


@pytest.mark.parametrize('element,tau', [('O', .5), ('Fe', .3)])
@pytest.mark.parametrize('above', [False, True])
def test_internal_and_boundary_counts_use_strict_original_weight_threshold(element, tau, above):
    weight = np.nextafter(tau, np.inf) if above else tau
    source, target = endpoint(('C', element), ((0, 1, weight),)), endpoint(('C', element), ())
    for partition in (((0, 1),), ((0,), (1,))):
        catalogue = saved(source, target, {0: 0, 1: 1},
                          pairs=tuple((atoms, atoms) for atoms in partition))
        decoded = decode_events(catalogue, event_slack=0)
        assert decoded.minimum_events == int(above)
        assert decoded.contributions(0)['total'] == int(above)


def test_partition_occurrences_remain_separate_and_full_catalogue_survives_export(ring_catalogue):
    catalogue = ring_catalogue
    catalogue.families = [catalogue.families[0]]
    next(iter(catalogue.branches.values())).clear()
    next(iter(catalogue.branches.values())).add(0)
    whole, singles = ((0, 1, 2, 3),), ((0,), (1,), (2,), (3,))
    catalogue.coupled_branches[0] = {tuple((atoms, atoms) for atoms in singles)}
    catalogue.families[0].provenance = [
        {'source_partition': whole, 'cuts': ()}, {'source_partition': singles, 'cuts': ()}]
    before = json.loads(json.dumps(catalogue.to_record()))
    decoded = decode_events(catalogue, event_slack=0)
    whole_id = decoded.fragment_partitions[0].index(whole)
    singles_id = decoded.fragment_partitions[0].index(singles)
    assert decoded.window.retained_partitions[0] == (whole_id,)
    assert decoded.attainable_counts(0, partition=whole_id) == (0,)
    with pytest.raises(ValueError, match='window'):
        decoded.attainable_counts(0, partition=singles_id)
    with pytest.raises(ValueError, match='window'):
        decoded.realize(0, partition=singles_id)
    record = json.loads(json.dumps(decoded_events_record(decoded)))
    assert record['catalogue'] == before
    assert record['schema'] == 'mappa.decoded_events/v10'
    assert record['model'] == 'independent_fragments_with_sweep'
    assert record['approximate'] is True and record['minimum_proven'] is False
    assert record['model_minimum_proven'] is True
    assert record['global_minimum_proven'] is False
    restored = FinalBranchCatalogue.from_record(catalogue.problem, record['catalogue'])
    assert json.loads(json.dumps(restored.to_record())) == before


def test_shared_fragment_tables_and_native_batches_match_correlated_scalar_counts(ring_catalogue):
    catalogue = ring_catalogue
    catalogue.families = [catalogue.families[0], replace(catalogue.families[0], provenance=[])]
    next(iter(catalogue.branches.values())).clear()
    next(iter(catalogue.branches.values())).update((0, 1))
    decoded = decode_events(catalogue, event_slack=None)
    assert len(decoded.inventory.fragments) == 1
    assert decoded.tables(0) is decoded.tables(1)
    choices = list(two_sided_choices(catalogue.problem.reactant.wbo))
    expected = [approximate_counts(catalogue.problem, dict(catalogue.families[0].mapping),
                ((0, 1, 2, 3),), mapping) for _, mapping in choices]
    for family in range(2):
        result = list(decoded.counts(family, (choice for choice, _ in choices), categories=True))
        assert [tuple(map(int, row)) for row in result] == [tuple(row.values()) for row in expected]


def test_original_wbo_is_used_without_changing_binary_matching_constraints():
    binary = endpoint(('V', 'O', 'O'), ((0, 1, 1.), (0, 2, 1.)))
    raw = AAMProblem(endpoint(('V', 'O', 'O'), ((0, 1, 2.2), (0, 2, .7))),
                     endpoint(('V', 'O', 'O'), ((0, 1, .8), (0, 2, 2.1))))
    identity = {0: 0, 1: 1, 2: 2}
    fragment = saved_fragment(AAMProblem(binary, binary), accepted_edges=((0, 1), (0, 2)))
    catalogue = saved(binary, binary, identity, fragments=(fragment,))
    catalogue.families.append(saved(binary, binary, identity).families[0])
    next(iter(catalogue.branches.values())).add(1)
    assert decode_events(catalogue, event_slack=0).retained_family_ids == (0, 1)
    decoded = decode_events(catalogue, event_problem=raw, event_slack=0)
    assert decoded.minimum_events == 0 and decoded.retained_family_ids == (0,)
    choice = ({'discovery': {0: 0, 1: 2, 2: 1}},)
    assert decoded.contributions(0, choice)['total'] == 0
    with pytest.raises(ValueError, match='window'):
        decoded.realize(0)
    record = json.loads(json.dumps(decoded_events_record(decoded)))
    assert record['matching_problem_sha256'] != record['event_problem_sha256']
    assert record['catalogue'] == json.loads(json.dumps(catalogue.to_record()))
    restored = FinalBranchCatalogue.from_record(problem_from_record(record['problem']), record['catalogue'])
    result = decode_events(restored, event_problem=problem_from_record(record['event_problem']), event_slack=0)
    assert result.attainable_counts(0) == (0,)


def test_nontransitive_local_relationship_is_not_replaced_by_independent_pair_choices():
    ep = endpoint('CCC', ((0, 1, 1.), (0, 2, 1.6), (1, 2, 2.2)))
    fragment = saved_fragment(AAMProblem(ep, ep), domains=tuple((a, (a,)) for a in range(3)),
        accepted_edges=((0, 1), (0, 2), (1, 2)), self_tolerance=.7)
    catalogue = saved(ep, ep, {0: 0, 1: 1, 2: 2}, fragments=(fragment,))
    allowed = self_permutations(ep.wbo, .7)
    assert len(allowed) == 3
    oracle = list(two_sided_choices(ep.wbo, tolerance=.7))
    bins = tuple(sorted({exact_count(catalogue.problem, mapping) for _, mapping in oracle}))
    decoded = decode_events(catalogue, event_slack=None)
    assert decoded.attainable_counts(0) == bins
    for choice, mapping in oracle:
        assert decoded.contributions(0, choice)['total'] == exact_count(catalogue.problem, mapping)
    forbidden = next(p for p in permutations(range(3)) if p not in allowed)
    with pytest.raises(ValueError):
        decoded.realize(0, ({'source_permutation': forbidden},))


def test_partial_metal_boundaries_are_fixed_and_wholly_unmapped_material_is_ignored():
    binary = endpoint(('Fe', 'Fe', 'C', 'C', 'C'),
                      tuple((a, b, 1.) for a, b in combinations(range(5), 2)))
    fragment = saved_fragment(AAMProblem(binary, binary), source_atoms=(0, 1),
                              accepted_edges=((0, 1),))
    catalogue = saved(binary, binary, {0: 0, 1: 1, 2: 2}, fragments=(fragment,),
                      pairs=(((0, 1), (0, 1)), ((2,), (2,))))
    raw = AAMProblem(
        endpoint(('Fe', 'Fe', 'C', 'C', 'C'), ((0, 2, .4), (1, 2, .9), (3, 4, 2.))),
        endpoint(('Fe', 'Fe', 'C', 'C', 'C'),
                 ((0, 2, .9), (1, 2, .4), (0, 3, .35), (1, 3, .7), (2, 3, .5), (3, 4, 3.))))
    decoded = decode_events(catalogue, event_problem=raw)
    assert decoded.minimum_events == 6 and decoded.attainable_counts(0) == (6,)
    choices = (None, ({'discovery': {0: 1, 1: 0}},))
    for choice in choices:
        assert decoded.contributions(0, choice)['counts'] == dict(
            breaking=2, forming=4, strengthening=0, weakening=0)
    assert [tuple(map(int, row)) for row in decoded.counts(0, choices, categories=True)] == [
        (2, 4, 0, 0), (2, 4, 0, 0)]
