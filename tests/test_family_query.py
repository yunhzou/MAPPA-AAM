from dataclasses import replace
from itertools import permutations, product
import numpy as np

from mappa import AAMProblem,AAMSearchConfig,search_aam,query_path
from mappa.mappa_core.domain import MolecularEndpoint
from mappa.mappa_core.core.matching.representation.records import Fragment
from mappa.mappa_core.core.symmetry import DiscoveryEdge, RelaxedIsomorphism, relaxed_automorphism
from mappa.mappa_core.post_processing.final_branches import FinalFamily


def endpoint(elements,edges):
    w=np.zeros((len(elements),len(elements)))
    for a,b in edges:w[a,b]=w[b,a]=1
    return MolecularEndpoint(tuple(elements),np.zeros((len(elements),3)),w)


def fragment_path(problem, *, atoms=None, targets=None, anchors=(), accepted=(),
                  self_tolerance=0.):
    atoms = tuple(range(problem.source_atom_count)) if atoms is None else atoms
    targets = atoms if targets is None else targets
    source = relaxed_automorphism(problem.reactant.elements, problem.reactant.wbo,
        atoms=atoms, anchors=tuple(r for r, _ in anchors), tolerance=self_tolerance)
    target = relaxed_automorphism(problem.product.elements, problem.product.wbo,
        atoms=targets, anchors=tuple(p for _, p in anchors), tolerance=self_tolerance)
    witness = tuple(zip(atoms, targets))
    domains = tuple((r, tuple(p for p in targets
        if problem.reactant.elements[r] == problem.product.elements[p])) for r in atoms)
    edges = tuple(DiscoveryEdge((a, b), problem.reactant.wbo[a, b], 0.) for a, b in accepted)
    target_edges = tuple((p, q) for p in targets for q in targets
                         if p < q and problem.product.wbo[p, q] >= .2)
    fragment = Fragment(RelaxedIsomorphism(source, target, witness, anchors,
        domains, edges, target_edges, .2))
    return FinalFamily(witness, accepted, (('relation', fragment),), (.2, 0.)).as_path(problem)


def test_open_pools_keep_same_pair_and_different_pair_assignments():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOOO',[(0,1),(2,3)]))
    result=search_aam(problem,AAMSearchConfig(seed_count=1))
    paths=list(result.graph.paths())
    # Compressed placements jointly cover both atoms in one pair and atoms
    # in different pairs, retaining their shared component choices.
    # Exhaustive oracle only for this tiny unit test, never production search.
    for image in permutations(range(4),2):
        reference=dict(enumerate(image))
        answers=[query_path(path,problem,reference,source_atoms=(0,1),timeout_ms=2000)
                 for path in paths]
        assert all(status in ('recovered','not_recovered') for status,_ in answers)
        witnesses=[witness for status,witness in answers if status=='recovered']
        assert len(witnesses)==1
        assert dict(witnesses[0]['mapping'])==reference


def test_correlated_fragment_is_not_independent_atom_membership():
    molecule=endpoint('CCCC',[(0,1),(1,2),(2,3)])
    problem=AAMProblem(molecule,molecule)
    result=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=.5))
    paths=list(result.graph.paths())
    for image in [(0,1,2,3),(3,2,1,0)]:
        assert any(query_path(path,problem,dict(enumerate(image)),source_atoms=tuple(range(4)),
            target_generators=((3,2,1,0),),timeout_ms=2000)[0]=='recovered' for path in paths)
    assert all(query_path(path,problem,{0:0,1:2,2:1,3:3},source_atoms=tuple(range(4)),
        target_generators=((3,2,1,0),),timeout_ms=2000)[0]=='not_recovered' for path in paths)


def test_full_hydrogen_mapping_is_returned_and_one_to_one():
    molecule=endpoint('OHOH',[(0,1),(2,3)])
    problem=AAMProblem(molecule,molecule)
    result=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10))
    # Check the complete compressed union without assuming which witness
    # supplies the swapped-component placement.
    answers=[query_path(path,problem,{0:2,2:0},source_atoms=(0,2),timeout_ms=3000)
             for path in result.graph.paths()]
    assert all(status in ('recovered','not_recovered') for status,_ in answers)
    witnesses=[witness for status,witness in answers if status=='recovered']
    assert witnesses
    witness=witnesses[0]
    mapping=dict(witness['mapping'])
    assert len(mapping)==len(set(mapping.values()))==4
    assert molecule.wbo[mapping[0],mapping[1]]==molecule.wbo[mapping[2],mapping[3]]==1


def test_one_target_cannot_supply_two_reference_atoms():
    molecule=endpoint('OO',[]);problem=AAMProblem(molecule,molecule)
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1)).graph.paths())
    assert query_path(path,problem,{0:0,1:0},source_atoms=(0,1),timeout_ms=2000)[0]=='not_recovered'


def test_fragment_choices_do_not_move_anchors():
    problem=AAMProblem(endpoint('OOO',[]),endpoint('OOO',[]))
    path=fragment_path(problem, anchors=((0, 0),))
    assert query_path(path,problem,{0:0,1:2},source_atoms=(0,1))[0]=='recovered'
    assert query_path(path,problem,{0:2,1:0},source_atoms=(0,1))[0]=='not_recovered'


def test_fixed_fragment_image_cannot_reach_an_extra_target_atom():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOO',[]))
    path=fragment_path(problem)
    assert query_path(path,problem,{0:0,1:2},source_atoms=(0,1))[0]=='not_recovered'


def test_placement_group_reaches_other_atom_set_without_mixing_its_atoms():
    problem = AAMProblem(endpoint('CC', [(0, 1)]),
                         endpoint('CCCC', [(0, 1), (2, 3)]))
    path = fragment_path(problem, atoms=(0, 1), targets=(0, 1), accepted=((0, 1),))
    edge = replace(path.graph.transitions[0], placement_generators=((2, 3, 0, 1),))
    path = replace(path, graph=replace(path.graph, transitions=(edge,)))
    expected = {(0, 1), (1, 0), (2, 3), (3, 2)}
    for images in permutations(range(4), 2):
        status, result = query_path(path, problem, dict(enumerate(images)), source_atoms=(0, 1))
        assert status == ('recovered' if images in expected else 'not_recovered')
        if status == 'recovered':
            assert dict(result['mapping']) == dict(enumerate(images))
            assert [action['kind'] for action in result['actions']] == ['relation', 'group']


def test_group_only_transition_preserves_its_joint_placement_choice():
    problem = AAMProblem(endpoint('CC', [(0, 1)]),
                         endpoint('CCCC', [(0, 1), (2, 3)]))
    path = fragment_path(problem, atoms=(0, 1), targets=(0, 1), accepted=((0, 1),))
    edge = replace(path.graph.transitions[0], match=None,
                   placement_generators=((2, 3, 0, 1),))
    path = replace(path, graph=replace(path.graph, transitions=(edge,)))
    reference = {0: 2, 1: 3}
    status, result = query_path(path, problem, reference, source_atoms=(0, 1))
    assert status == 'recovered' and dict(result['mapping']) == reference
    assert [action['kind'] for action in result['actions']] == ['group']
    assert result['checked_fragment_edges'] == 0
    assert query_path(path, problem, {0: 3, 1: 2}, source_atoms=(0, 1))[0] == 'not_recovered'


def test_subgraph_placement_replay_extends_by_identity_in_enclosing_problem():
    from mappa import decode_events
    from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue

    problem = AAMProblem(endpoint('CCC', [(0, 1), (1, 2)]),
                         endpoint('CCCO', [(0, 1), (1, 2)]))
    path = fragment_path(problem, atoms=(0, 1, 2), targets=(0, 1, 2),
                         anchors=((0, 0), (1, 1), (2, 2)), accepted=((0, 1), (1, 2)))
    edge = replace(path.graph.transitions[0], placement_generators=((2, 1, 0),))
    context = replace(path.context, target_atoms=(0, 1, 2))
    path = replace(path, graph=replace(path.graph, transitions=(edge,), contexts=(context,)))
    reference = {0: 2, 1: 1, 2: 0}
    status, witness = query_path(path, problem, reference, source_atoms=(0, 1, 2))
    assert status == 'recovered' and dict(witness['mapping']) == reference
    assert dict(witness['actions'][-1]['permutation'])[3] == 3

    catalogue = FinalBranchCatalogue(problem).add_graph(path.graph)
    family = catalogue.families[0]
    assert family.actions[-1] == ('group', ((2, 1, 0, 3),))
    decoded = decode_events(catalogue)
    assert decoded.realize(0, (None, (1,))) == reference
    assert dict(family.as_path(problem).realize(placement_choices={0: (1,)}).mapping) == reference


def test_placement_program_keeps_later_choices_correlated_with_earlier_moves():
    from mappa.mappa_core.core.search.graph import (
        AAMSearchGraph, FragmentTransition, SearchContext, SearchState, SearchStop,
    )
    from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments
    source = endpoint('CCCC', [(0, 1), (2, 3)])
    target = endpoint('CCCCCC', [(0, 1), (2, 3), (4, 5)])
    problem = AAMProblem(source, target)
    first = fragment_path(problem, atoms=(0, 1), targets=(0, 1),
                          accepted=((0, 1),)).fragments[0]
    second = fragment_path(problem, atoms=(2, 3), targets=(2, 3),
                           accepted=((2, 3),)).fragments[0]
    rotation = (2, 3, 4, 5, 0, 1)
    exchange = (0, 1, 4, 5, 2, 3)  # The second placement fixes earlier locks.
    graph = AAMSearchGraph(
        (SearchContext(tuple(range(4)), tuple(range(6)), (0, 2)),), (0,),
        (SearchState(0, 0, (), (), ()),
         SearchState(1, 0, ((0, 0), (1, 1)), ((0, 0), (1, 0)), ()),
         SearchState(2, 0, tuple(enumerate(range(4))),
                     ((0, 0), (1, 0), (2, 1), (3, 1)), ())),
        (FragmentTransition(0, 0, 1, 0, (0,), first, 0, (rotation,)),
         FragmentTransition(1, 1, 2, 2, (1,), second, 1, (exchange,))),
        (SearchStop(2, 'objective_met'),))
    path = next(graph.paths())
    # Independent tiny oracle: orient each pair, then place the two distinct
    # pairs into any two of the three components. There are 6 * 2 * 2 choices.
    expected = {
        left + right
        for a, b in permutations(range(3), 2)
        for left, right in product(permutations((2*a, 2*a+1)),
                                   permutations((2*b, 2*b+1)))
    }
    projected = project_path_assignments(path, problem, range(4), assignment_limit=24)
    assert {tuple(mapping[r] for r in range(4)) for mapping in projected} == expected
    reference = {0: 5, 1: 4, 2: 0, 3: 1}
    status, witness = query_path(path, problem, reference, source_atoms=tuple(range(4)))
    assert status == 'recovered' and dict(witness['mapping']) == reference
    assert [action['kind'] for action in witness['actions']] == ['relation', 'group'] * 2
    actual = path.mapping
    for action in witness['actions']:
        transport = dict(action['permutation'])
        actual = {r: transport.get(p, p) for r, p in actual.items()}
    assert actual == reference
    assert query_path(path, problem, {0: 0, 1: 2, 2: 1, 3: 3},
                      source_atoms=tuple(range(4)))[0] == 'not_recovered'


def test_reference_heavy_assignment_still_requires_feasible_hydrogens():
    problem=AAMProblem(endpoint('OOH',[(0,2)]),endpoint('OOH',[(0,2)]))
    path=fragment_path(problem, accepted=((0, 2),))
    assert query_path(path,problem,{0:0},source_atoms=(0,))[0]=='recovered'
    assert query_path(path,problem,{0:1},source_atoms=(0,))[0]=='not_recovered'


def test_reference_propagation_keeps_correlated_permutation_choices():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
    for target in permutations(range(3)):
        constraints=FiniteConstraints()
        values=constraints.act([constraints.constant(i) for i in range(3)],[(1,2,0)],'ring')
        for value,wanted in zip(values,target):constraints.restrict(value[0],{wanted})
        solver=z3.Solver();constraints.emit(solver)
        expected=target in ((0,1,2),(1,2,0),(2,0,1))
        assert (solver.check()==z3.sat)==expected
        assert (constraints.solve(2000)[0]==z3.sat)==expected


def test_reference_propagation_dynamic_lookup_preserves_all_solutions():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
    # The reference fixes the output, not a freely chosen input slot.
    for selector in (0,1):
        for left in (2,3):
            for right in (2,3):
                constraints=FiniteConstraints()
                choice=constraints.symbol((0,1))
                a=constraints.symbol((2,3));b=constraints.symbol((2,3))
                output=constraints.lookup(choice,{0:a,1:b})
                constraints.restrict(output[0],{2})
                constraints.restrict(choice[0],{selector})
                constraints.restrict(a[0],{left});constraints.restrict(b[0],{right})
                solver=z3.Solver();constraints.emit(solver)
                assert (solver.check()==z3.sat)==((left,right)[selector]==2)


def test_reference_propagation_overlapping_pools_matches_exhaustive_oracle():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
    # Two correlated permutations on overlapping target pools, with one slot
    # fixed in the second operation. Compare every output, including negatives.
    expected=set()
    for first in permutations((0,1,2)):
        for second in permutations((1,2,3)):
            table=dict(zip((1,2,3),second))
            if table[3]!=3:continue
            expected.add(tuple(table.get(v,v) for v in (*first,3)))
    for target in permutations(range(4)):
        constraints=FiniteConstraints()
        values=[constraints.constant(i) for i in range(4)]
        for pool,fixed in (((0,1,2),()),((1,2,3),(3,))):
            table={p:constraints.symbol(pool) for p in pool}
            constraints.distinct(v[0] for v in table.values())
            for p in fixed:constraints.restrict(table[p][0],{p})
            values=[constraints.lookup(v,table) for v in values]
        for value,wanted in zip(values,target):constraints.restrict(value[0],{wanted})
        solver=z3.Solver();constraints.emit(solver)
        assert (solver.check()==z3.sat)==(target in expected)
        assert (constraints.solve(2000)[0]==z3.sat)==(target in expected)


def test_reference_atom_set_rejects_before_action_construction(monkeypatch):
    problem=AAMProblem(endpoint('CCC',[]),endpoint('CCC',[]))
    path=fragment_path(problem, atoms=(0, 1))
    def forbidden(*args):raise AssertionError('Atom-set mismatch needs no group query')
    monkeypatch.setattr('mappa.mappa_core.core.symmetry.constraints.group_factors',forbidden)
    assert query_path(path,problem,{0:0,1:2},source_atoms=(0,1))[0]=='not_recovered'


def test_reference_rejects_cross_fragment_images_before_relation_compilation(monkeypatch):
    molecule = endpoint('CCCC', [(0, 1), (2, 3)])
    problem = AAMProblem(molecule, molecule)
    fragments = tuple(fragment_path(problem, atoms=pair, accepted=(pair,)).fragments[0]
                      for pair in ((0, 1), (2, 3)))
    path = FinalFamily(tuple(enumerate(range(4))), ((0, 1), (2, 3)),
        tuple(('relation', fragment) for fragment in fragments), (.2, 0.)).as_path(problem)
    options = dict(source_atoms=tuple(range(4)),
                   source_generators=((1, 0, 3, 2),),
                   target_generators=((1, 0, 3, 2),))
    assert query_path(path, problem, {0: 1, 1: 0, 2: 3, 3: 2}, **options)[0] == 'recovered'

    def forbidden(*args, **kwargs):
        raise AssertionError('Separate fragment images already exclude this reference')

    monkeypatch.setattr('mappa.mappa_core.core.symmetry.constraints.compile_relaxed_relation',
                        forbidden)
    # All requested targets occur in the family, but neither local rearrangement
    # nor the declared endpoint equivalence can move an atom between these pairs.
    assert query_path(path, problem, {0: 2, 1: 3, 2: 0, 3: 1}, **options)[0] == 'not_recovered'


def test_reference_endpoint_equivalence_preserves_partial_and_unmatched_atoms(monkeypatch):
    from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments

    molecule = endpoint('CCCC', [(0, 1), (2, 3)])
    problem = AAMProblem(molecule, molecule)
    path = fragment_path(problem, atoms=(0, 1), accepted=((0, 1),))
    exchange = ((2, 3, 0, 1),)
    options = dict(source_atoms=tuple(range(4)), source_generators=exchange,
                   target_generators=exchange)
    reference = {2: 2, 3: 3}
    status, witness = query_path(path, problem, reference, **options)
    assert status == 'recovered'
    physical = dict(witness['mapping'])
    assert set(physical) == {0, 1} and set(physical.values()) == {0, 1}
    for missing in ('source_generators', 'target_generators'):
        assert query_path(path, problem, reference,
                          **{key: value for key, value in options.items() if key != missing})[0] == 'not_recovered'
    assert query_path(path, problem, {2: 2}, complete_reference=False, **options)[0] == 'recovered'
    assert query_path(path, problem, {2: 2}, **options)[0] == 'not_recovered'
    assert query_path(path, problem, {0: 0, 1: 1, 2: 4, 3: 4},
                      source_atoms=tuple(range(4)))[0] == 'recovered'
    assert project_path_assignments(path, problem, (0, 1), assignment_limit=2) == [
        {0: 0, 1: 1}, {0: 1, 1: 0}]
    assert project_path_assignments(path, problem, (0, 2), assignment_limit=2) == []

    def forbidden(*args, **kwargs):
        raise AssertionError('A required unmatched atom is already known to be mapped')

    monkeypatch.setattr('mappa.mappa_core.core.symmetry.constraints.compile_relaxed_relation',
                        forbidden)
    assert query_path(path, problem, {0: 0}, source_atoms=tuple(range(4)))[0] == 'not_recovered'


def test_reference_bounds_include_both_self_actions_around_fixed_discovery():
    from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments

    problem = AAMProblem(endpoint('CCC', [(0, 1), (0, 2)]),
                         endpoint('CCC', [(0, 1), (1, 2)]))
    path = fragment_path(problem, accepted=((0, 1),))
    fragment = Fragment(replace(path.fragments[0].relation,
                                domains=tuple((atom, (atom,)) for atom in range(3))))
    path = FinalFamily(tuple(enumerate(range(3))), ((0, 1),),
                       (('relation', fragment),), (.2, 0.)).as_path(problem)
    # Source leaves 1/2 may exchange; target leaves 0/2 may exchange. Either
    # exchange alone misses some of the four physically permitted mappings.
    expected = {(0, 1, 2), (0, 2, 1), (2, 1, 0), (2, 0, 1)}
    for images in permutations(range(3)):
        status, witness = query_path(path, problem, dict(enumerate(images)),
                                     source_atoms=(0, 1, 2))
        assert status == ('recovered' if images in expected else 'not_recovered')
        if status == 'recovered':
            assert dict(witness['mapping']) == dict(enumerate(images))
    projected = project_path_assignments(path, problem, range(3), assignment_limit=4)
    assert {tuple(mapping[r] for r in range(3)) for mapping in projected} == expected


def test_reference_and_projection_keep_nontransitive_two_sided_relations():
    from mappa.mappa_workflow.mapping_analytics.family import project_path_assignments

    molecule = endpoint('CCCC', [(0, 1), (0, 2), (0, 3)])
    weights = molecule.wbo.copy()
    for atom, weight in enumerate((1., 1.6, 2.2), 1):
        weights[0, atom] = weights[atom, 0] = weight
    molecule = MolecularEndpoint(molecule.elements, molecule.coordinates, weights)
    problem = AAMProblem(molecule, molecule)
    path = fragment_path(problem, accepted=((0, 1), (0, 2), (0, 3)), self_tolerance=.7)
    fragment = Fragment(replace(path.fragments[0].relation,
                                domains=tuple((atom, (atom,)) for atom in range(4))))
    path = FinalFamily(tuple(enumerate(range(4))), ((0, 1), (0, 2), (0, 3)),
                       (('relation', fragment),), (.2, 0.)).as_path(problem)
    # Exact discovery fixes each weighted arm. Each self-side may swap adjacent
    # weights, once. Together they give five permutations, not their six-element closure.
    expected = {(0, *tail) for tail in permutations((1, 2, 3)) if tail != (3, 2, 1)}
    for images in permutations(range(4)):
        status, witness = query_path(path, problem, dict(enumerate(images)),
                                     source_atoms=tuple(range(4)))
        assert status == ('recovered' if images in expected else 'not_recovered')
        if status == 'recovered':
            assert dict(witness['mapping']) == dict(enumerate(images))
    projected = project_path_assignments(path, problem, range(4), assignment_limit=5)
    assert {tuple(mapping[r] for r in range(4)) for mapping in projected} == expected

    # Reference equivalence is chemical endpoint normalization, not a further
    # physical self rearrangement. Either endpoint may normalize the sixth reference.
    reference = {0: 0, 1: 3, 2: 2, 3: 1}
    for side in ('source_generators', 'target_generators'):
        status, witness = query_path(path, problem, reference, source_atoms=tuple(range(4)),
                                     **{side: ((0, 3, 2, 1),)})
        assert status == 'recovered'
        assert tuple(dict(witness['mapping'])[r] for r in range(4)) in expected


def test_reference_joint_rejection_does_not_build_output_expressions(monkeypatch):
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints

    problem = AAMProblem(endpoint('CC', [(0, 1)]),
                         endpoint('CCCC', [(0, 1), (2, 3)]))
    path = fragment_path(problem, accepted=((0, 1),))
    edge = replace(path.graph.transitions[0], match=None,
                   placement_generators=((2, 3, 0, 1),))
    path = replace(path, graph=replace(path.graph, transitions=(edge,)))
    expression = FiniteConstraints.expression

    def require_live_query(self, value):
        assert not self.failed, 'A rejected reference needs no output expressions'
        return expression(self, value)

    monkeypatch.setattr(FiniteConstraints, 'expression', require_live_query)
    # Each requested image is individually reachable. Their shared component
    # choice cannot stay for atom 0 and move for atom 1 at the same time.
    assert query_path(path, problem, {0: 0, 1: 3},
                      source_atoms=(0, 1))[0] == 'not_recovered'


def test_reference_respects_complete_fragment_ring_correlation():
    edges=((0, 1), (1, 2), (2, 3), (0, 3))
    problem=AAMProblem(endpoint('CCCC',edges),endpoint('CCCC',edges))
    path=fragment_path(problem, accepted=edges)
    for expected,reference in [('recovered',{0:1,1:2,2:3,3:0}),
                               ('not_recovered',{0:1,1:0,2:2,3:3})]:
        status,result=query_path(path,problem,reference,source_atoms=tuple(range(4)))
        assert status==expected
        if status=='recovered':
            assert dict(result['mapping'])==reference


def test_independent_fragment_constraint_cache_and_joint_conflict():
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, _component_cache
    _component_cache.clear()
    constraints=FiniteConstraints()
    for _ in range(2):
        values=[constraints.symbol(range(3))[0] for _ in range(3)]
        for a,b in ((0,1),(1,2),(2,0)):
            constraints.relation((values[a],values[b]),tuple((x,y) for x in range(3) for y in range(3) if x!=y))
    status,_,stats=constraints.solve(2000)
    assert status==z3.sat
    assert stats['components']==2 and stats['cache_hits']==1 and stats['solver_calls']==1
    # Sharing a choice means these constraints must not be accepted separately.
    joint=FiniteConstraints();a,b,c=[joint.symbol((0,1))[0] for _ in range(3)]
    joint.relation((a,b),((0,0),(1,1)))
    joint.relation((b,c),((0,0),(1,1)))
    joint.relation((a,c),((0,1),(1,0)))
    assert joint.solve(2000)[0]==z3.unsat


def test_constraint_checkpoints_relocate_to_the_shared_core():
    import io
    import pickle
    import z3
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
    from mappa.mappa_core.utils.io.checkpoints import _load_checkpoint
    system = FiniteConstraints()
    first, second = system.symbol((0, 1)), system.symbol((0, 1))
    system.distinct((first[0], second[0]))
    saved = pickle.dumps(system, protocol=0).replace(
        b'cmappa.mappa_core.core.symmetry.constraints\nFiniteConstraints\n',
        b'cmappa.mappa_workflow.mapping_analytics.constraints\nReferenceConstraints\n').replace(
        b'cmappa.mappa_core.core.symmetry.constraints\nChoice\n',
        b'cmappa.mappa_workflow.mapping_analytics.constraints\nChoice\n')
    restored = _load_checkpoint(io.BytesIO(saved))
    assert isinstance(restored, FiniteConstraints)
    assert restored.solve(2000)[0] == z3.sat


def test_tolerant_symmetry_retains_all_small_query_mappings():
    mol=endpoint('CCCHHH',[(0,1),(1,2),(0,2),(0,3),(1,4),(2,5)])
    weights=mol.wbo.copy()
    weights[0,1]=weights[1,0]=2
    mol=MolecularEndpoint(mol.elements,mol.coordinates,weights)
    problem=AAMProblem(mol,mol)
    aam=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=1))
    paths=list(aam.graph.paths())
    # Every tolerance-compatible placement must remain in the union.
    from mappa import query_path
    for p in permutations(range(3)):
        mapping=dict(enumerate((*p,*(i+3 for i in p))))
        assert any(query_path(path,problem,mapping,source_atoms=tuple(range(6)),
                              timeout_ms=2000)[0]=='recovered' for path in paths)
