"""Independent membership and replay oracles for shared fragment verification."""
from itertools import permutations
from dataclasses import replace
from types import SimpleNamespace

from mappa import AAMProblem, AAMSearchConfig, search_aam, query_path
from mappa.mappa_core.core.search.queries.reference import ReferenceQuery
from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, FinalFamily
from mappa.mappa_workflow.mapping_analytics.queries import query_catalogue
from test_family_query import endpoint, fragment_path


def test_shared_fragments_are_checked_once_then_composed_and_replayed():
    molecule = endpoint('CCCCC', [(0, 1), (2, 3), (3, 4)])
    problem = AAMProblem(molecule, molecule)
    first = fragment_path(problem, atoms=(0, 1), accepted=((0, 1),)).fragments[0]
    second = fragment_path(problem, atoms=(2, 3, 4), accepted=((2, 3), (3, 4))).fragments[0]
    actions = (('relation', first), ('relation', second))
    reference = {0: 1, 1: 0, 2: 4, 3: 3, 4: 2}
    query = ReferenceQuery(problem, reference)
    for program in (actions, actions[::-1]):
        status, evidence = query.query(tuple(enumerate(range(5))), program)
        assert status == 'recovered'
        assert dict(evidence['mapping']) == reference
    assert query.stats['fragment_queries'] == 2
    assert query.stats['cache_hits'] == 2
    assert query.stats['joint_queries'] == 0
    family = FinalFamily(tuple(enumerate(range(5))), (), actions, (.2, 1.))
    result = query_catalogue(SimpleNamespace(problem=problem, families=[family]), reference)
    proof = result['evidence']
    replay = family.as_path(problem).realize(fragment_choices=proof['fragment_choices'],
                                            placement_choices=proof['placement_choices'])
    assert dict(replay.mapping) == reference


def test_impossible_composition_needs_no_fragment_solver():
    molecule = endpoint('CCCC', [(0, 1), (2, 3)])
    problem = AAMProblem(molecule, molecule)
    actions = tuple(('relation', fragment_path(problem, atoms=pair,
        accepted=(pair,)).fragments[0]) for pair in ((0, 1), (2, 3)))
    query = ReferenceQuery(problem, {0: 2, 2: 0})
    assert query.query(tuple(enumerate(range(4))), actions)[0] == 'not_recovered'
    assert query.stats['composition_rejections'] == 1
    assert query.stats['fragment_queries'] == query.stats['joint_queries'] == 0


def test_coupled_placements_do_not_accept_individually_reachable_atoms():
    molecule = endpoint('CCC', [(0, 1), (1, 2), (2, 0)])
    problem = AAMProblem(molecule, molecule)
    family = FinalFamily(tuple(enumerate(range(3))), (),
                         (('group', ((1, 2, 0),)),), (.2, 1.))
    catalogue = SimpleNamespace(problem=problem, families=[family])
    for image in permutations(range(3)):
        reference = dict(enumerate(image))
        result = query_catalogue(catalogue, reference)
        assert result['status'] == ('recovered' if image in
            ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else 'not_recovered')
        if result['status'] == 'recovered':
            assert dict(result['evidence']['mapping']) == reference


def test_heavy_reference_keeps_hydrogen_feasibility_and_anchors():
    molecule = endpoint('OOH', [(0, 2)])
    problem = AAMProblem(molecule, molecule)
    fragment = fragment_path(problem, accepted=((0, 2),)).fragments[0]
    family = FinalFamily(fragment.witness, (), (('relation', fragment),), (.2, 1.))
    catalogue = SimpleNamespace(problem=problem, families=[family])
    assert query_catalogue(catalogue, {0: 1})['status'] == 'not_recovered'
    assert query_catalogue(catalogue, {0: 0})['status'] == 'recovered'
    anchored = fragment_path(problem, anchors=((0, 0),)).fragments[0]
    family.actions = (('relation', anchored),)
    assert query_catalogue(catalogue, {0: 1})['status'] == 'not_recovered'


def test_unannotated_hydrogen_placement_does_not_join_heavy_fragments():
    molecule = endpoint('OOHHHH', [(0, 2), (0, 3), (1, 4), (1, 5)])
    problem = AAMProblem(molecule, molecule)
    fragment = fragment_path(problem, atoms=(0, 1), accepted=()).fragments[0]
    actions = (('relation', fragment), ('group', ((0, 1, 4, 5, 2, 3),)))
    query = ReferenceQuery(problem, {0: 1, 1: 0})
    status, proof = query.query(tuple(enumerate(range(6))), actions)
    assert status == 'recovered' and dict(proof['mapping']) == {0: 1, 1: 0, 2: 2, 3: 3, 4: 4, 5: 5}
    assert query.stats['joint_queries'] == 0


def test_small_catalogues_equal_exact_saved_program_union():
    # Exhaust every reference on independent small structures. This includes
    # overlapping local/placement actions, unequal endpoint sizes and hydrogens.
    problems = [
        AAMProblem(endpoint('CC', []), endpoint('CCCC', [(0, 1), (2, 3)])),
        AAMProblem(endpoint('CCCC', [(0, 1), (1, 2), (2, 3)]),
                   endpoint('CCCC', [(0, 1), (1, 2), (2, 3), (3, 0)])),
        AAMProblem(endpoint('OHOH', [(0, 1), (2, 3)]),
                   endpoint('OHOH', [(0, 1), (2, 3)])),
    ]
    for problem in problems:
        result = search_aam(problem, AAMSearchConfig(seed_count=1, cut_floor=10))
        catalogue = FinalBranchCatalogue(problem).add_graph(result.graph)
        atoms = tuple(r for r, element in enumerate(problem.reactant.elements) if element != 'H')
        for images in permutations(range(problem.target_atom_count), len(atoms)):
            reference = dict(zip(atoms, images))
            answers = [query_path(family.as_path(problem), problem, reference,
                source_atoms=atoms, timeout_ms=5000)[0] for family in catalogue.families]
            assert 'unknown' not in answers
            expected = 'recovered' if 'recovered' in answers else 'not_recovered'
            answer = query_catalogue(catalogue, reference, timeout_ms=5000)
            assert answer['status'] == expected
            if expected == 'recovered':
                actual = dict(answer['evidence']['mapping'])
                assert all(actual[r] == p for r, p in reference.items())
                assert len(set(actual.values())) == len(actual)


def test_both_self_actions_remain_allowed_when_discovery_is_fixed():
    problem = AAMProblem(endpoint('CCC', [(0, 1), (0, 2)]),
                         endpoint('CCC', [(0, 1), (1, 2)]))
    original = fragment_path(problem, accepted=((0, 1),)).fragments[0]
    fragment = Fragment(replace(original.relation,
        domains=tuple((atom, (atom,)) for atom in range(3))))
    family = FinalFamily(fragment.witness, (), (('relation', fragment),), (.2, 0.))
    catalogue = SimpleNamespace(problem=problem, families=[family])
    expected = {(0, 1, 2), (0, 2, 1), (2, 1, 0), (2, 0, 1)}
    for images in permutations(range(3)):
        answer = query_catalogue(catalogue, dict(enumerate(images)))
        assert answer['status'] == ('recovered' if images in expected else 'not_recovered')


def test_overlapping_actions_keep_order_instead_of_group_closure():
    molecule = endpoint('CCC', [])
    problem = AAMProblem(molecule, molecule)
    actions = (('group', ((1, 0, 2),)), ('group', ((0, 2, 1),)))
    family = FinalFamily(tuple(enumerate(range(3))), (), actions, (.2, 1.))
    catalogue = SimpleNamespace(problem=problem, families=[family])
    expected = {(0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 0, 1)}
    for images in permutations(range(3)):
        answer = query_catalogue(catalogue, dict(enumerate(images)))
        assert answer['status'] == ('recovered' if images in expected else 'not_recovered')


def test_inconclusive_fragment_stays_unknown_and_can_be_reused(monkeypatch):
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
    molecule = endpoint('CCC', [])
    problem = AAMProblem(molecule, molecule)
    query = ReferenceQuery(problem, {0: 1})
    monkeypatch.setattr(FiniteConstraints, 'solve', lambda *args: (z3.unknown, None, 0))
    program = (('group', ((1, 2, 0),)),)
    for _ in range(2):
        assert query.query(tuple(enumerate(range(3))), program)[0] == 'unknown'
    assert query.stats['joint_queries'] == 1
    assert query.stats['cache_hits'] == 1


def test_random_ordered_compositions_against_explicit_small_permutations():
    import random
    rng = random.Random(42)
    molecule = endpoint('CCCC', [])
    problem = AAMProblem(molecule, molecule)
    identity = tuple(range(4))
    for _ in range(100):
        actions, expected = [], {identity}
        for _ in range(4):
            atoms = tuple(sorted(rng.sample(range(4), rng.choice((2, 3)))))
            if rng.randrange(2):
                fragment = fragment_path(problem, atoms=atoms).fragments[0]
                actions.append(('relation', fragment))
                choices = []
                for images in permutations(atoms):
                    permutation = list(identity)
                    for a, b in zip(atoms, images):
                        permutation[a] = b
                    choices.append(tuple(permutation))
            else:
                generator = list(identity)
                for a, b in zip(atoms, atoms[1:] + atoms[:1]):
                    generator[a] = b
                actions.append(('group', (tuple(generator),)))
                choices, current = [], identity
                for _ in atoms:
                    choices.append(current)
                    current = tuple(generator[p] for p in current)
            expected = {tuple(choice[p] for p in previous)
                        for previous in expected for choice in choices}
        for images in permutations(range(4)):
            reference = dict(enumerate(images))
            status, proof = ReferenceQuery(problem, reference).query(
                tuple(enumerate(identity)), actions)
            assert status == ('recovered' if images in expected else 'not_recovered')
            if status == 'recovered':
                assert dict(proof['mapping']) == reference
