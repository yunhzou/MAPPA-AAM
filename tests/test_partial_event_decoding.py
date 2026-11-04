"""One event rule for partial and full maps; exhaustive tiny graphs are the oracle."""
from mappa.mappa_workflow.mapping_analytics.queries import query_assignments
from itertools import combinations, permutations
import random

import numpy as np
import pytest

from mappa import AAMProblem, MolecularEndpoint
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_core.core.symmetry.event_counts import _event_model
from mappa.mappa_workflow.mapping_analytics.event_query import _same_event_orbit
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, FinalFamily
from mappa.mappa_core.post_processing.hierarchical_mapping import RefinementConfig
from mappa.mappa_core.post_processing.hierarchical_mapping import refine_mapping
from mappa.mappa_core.core.symmetry.family import compile_path
from reference_event_decoder import extract_path_events
from growth_cases import saved_fragment, represented


def endpoint(elements, edges):
    matrix = np.zeros((len(elements), len(elements)))
    for a, b, w in edges:
        matrix[a, b] = matrix[b, a] = w
    return MolecularEndpoint(tuple(elements), np.zeros((len(elements), 3)), matrix)


def oracle(problem, mapping, threshold=.5):
    inverse = {b: a for a, b in mapping.items()}
    result = []
    for source, target, correspondences in [(problem.reactant, problem.product, mapping),
                                             (problem.product, problem.reactant, inverse)]:
        count = 0
        for a, b in combinations(range(source.atom_count), 2):
            if a not in correspondences and b not in correspondences:
                continue
            other = target.wbo[correspondences[a], correspondences[b]] if a in correspondences and b in correspondences else 0
            count += source.wbo[a, b] - other > threshold
        result.append(int(count))
    return tuple(result)


def catalogue(problem, mapping, targets=None):
    """Represent each selected target set by its own completed relation."""
    cat = FinalBranchCatalogue(problem)
    source_atoms = tuple(sorted(mapping))
    target_sets = (tuple(sorted(mapping.values())),) if targets is None else combinations(targets, len(mapping))
    for target_atoms in target_sets:
        witness = mapping if targets is None else dict(zip(source_atoms, target_atoms))
        fragments = () if not mapping else (saved_fragment(problem,
            source_atoms=source_atoms, target_atoms=target_atoms, witness=witness,
            anchors=tuple(mapping.items()) if targets is None else ()),)
        cat.families.append(FinalFamily(tuple(sorted(witness.items())), (),
            tuple(('relation', fragment) for fragment in fragments), (.2, 1.)))
        pairs = ((source_atoms, tuple(target_atoms)),) if mapping else ()
        cat.branches[((.2, 1.), pairs)] = {len(cat.families) - 1}
    return cat


def test_disconnected_unmatched_material_is_ignored_in_both_directions():
    core = endpoint('C', [])
    with_water = endpoint('COHH', [(1, 2, 1), (1, 3, 1)])
    for r, p in [(core, with_water), (with_water, core)]:
        problem = AAMProblem(r, p)
        decoded = refine_mapping(catalogue(problem, {0: 0}), RefinementConfig(event_slack=0))
        assert decoded.complete and decoded.minimum_proven and decoded.observed_minimum == 0
        assert decoded.candidates[0].mapping == {0: 0}
        assert decoded.candidates[0].endpoint_events == {'broken': (), 'formed': ()}


def test_mapped_boundary_and_interfragment_bonds_use_same_rule():
    r = endpoint('CCOH', [(0,1,1), (1,2,1), (2,3,1)])
    p = endpoint('CCNH', [(0,2,1), (2,3,1)])
    problem = AAMProblem(r, p)
    decoded = refine_mapping(catalogue(problem, {0:0,1:1}), RefinementConfig(event_slack=0))
    assert decoded.minimum_proven and decoded.observed_minimum == 3
    candidate = decoded.candidates[0]
    assert candidate.endpoint_events == {'broken': ((0,1),(1,2)), 'formed': ((0,2),)}
    reverse = refine_mapping(catalogue(AAMProblem(p,r), {0:0,1:1}), RefinementConfig(event_slack=0))
    assert reverse.minimum_proven and reverse.observed_minimum == 3
    assert reverse.candidates[0].endpoint_events == {'broken': ((0,2),), 'formed': ((0,1),(1,2))}


@pytest.mark.parametrize('nr,np_,mapped', [(2,4,(0,1)), (4,3,(0,2)), (3,3,(0,2))])
def test_partial_symbolic_events_and_classes_match_all_injections(nr, np_, mapped):
    import time
    import z3
    rng = random.Random(71+nr+np_)
    for trial in range(3):
        def molecule(n):
            return endpoint('C'*n, [(a,b,rng.choice((0.,.4,.6,1.,1.8))) for a,b in combinations(range(n),2)])
        problem = AAMProblem(molecule(nr),molecule(np_))
        maps = [dict(zip(mapped, images)) for images in permutations(range(np_),len(mapped))]
        index = SignedEventIndex(problem)
        patterns = [index.describe(index.vector(m)) for m in maps]
        assert [tuple(index.counts([index.vector(m)])[0]) for m in maps] == [oracle(problem,m) for m in maps]
        cat = catalogue(problem, maps[0], range(np_))
        observed = set()
        for family in cat.families:
            path = family.as_path(problem)
            compiled = compile_path(path,problem,{},source_atoms=(),complete_reference=False)
            terms,total,_ = _event_model(compiled,index,time.perf_counter()+30)
            equalities = [_same_event_orbit(terms,total,p,index) for p in patterns[:3]]
            target_set = set(dict(family.mapping).values())
            for mapping, pattern in zip(maps,patterns):
                compiled.solver.push()
                compiled.solver.add(*(compiled.values[a][0]==p for a,p in mapping.items()))
                if set(mapping.values()) != target_set:
                    assert compiled.solver.check() == z3.unsat
                    compiled.solver.pop()
                    continue
                assert compiled.solver.check()==z3.sat
                assert compiled.solver.model().eval(total).as_long()==sum(oracle(problem,mapping))
                for expected,equality in zip(patterns[:3],equalities):
                    compiled.solver.push()
                    compiled.solver.add(equality)
                    assert (compiled.solver.check()==z3.sat)==(expected['id']==pattern['id'])
                    compiled.solver.pop()
                compiled.solver.pop()
            result = extract_path_events(path, problem, index, seconds=30, max_patterns=None)
            assert result['complete']
            observed.update(p['id'] for p in result['patterns'])
        assert observed == {p['id'] for p in patterns}


def test_partial_target_symmetry_needs_no_solver(monkeypatch):
    r = endpoint('C', [])
    p = endpoint('CCC', [(0,1,1),(1,2,1),(0,2,1)])
    cat = catalogue(AAMProblem(r,p),{0:0},range(3))
    def forbidden(*args,**kwargs):
        raise AssertionError('Exact target symmetry needs no solver')
    monkeypatch.setattr('mappa.mappa_core.core.symmetry.family.compile_path',forbidden)
    result=refine_mapping(cat,RefinementConfig(event_slack=0))
    assert result.minimum_proven and result.observed_minimum==2
    assert len(result.candidates)==1
    assert all(report['solver_queries'] == 0 for report in result.family_reports)


def test_saved_partial_graph_survives_catalogue_and_artifact():
    from mappa.mappa_core.utils.io.records import refined_mapping_record
    problem=AAMProblem(endpoint('CO',[(0,1,1)]),endpoint('CON',[(0,2,1)]))
    family=catalogue(problem,{0:0,1:1}).families[0]
    cat=FinalBranchCatalogue(problem).add_graph(family.as_path(problem).graph)
    assert len(cat.families)==1 and cat.families[0].mapping==family.mapping
    decoded=refine_mapping(cat,RefinementConfig(event_slack=0))
    record=refined_mapping_record(decoded)
    assert record['complete'] and record['observed_minimum']==2
    assert record['candidates'][0]['endpoint_events']['formed']==((0,2),)
    candidate=decoded.candidates[0]
    assert query_assignments(decoded, candidate,{0:0,1:1}).status=='allowed'


def test_impossible_assignment_needs_no_event_composition(monkeypatch):
    problem = AAMProblem(endpoint('CC', []), endpoint('CCC', []))
    decoded = refine_mapping(catalogue(problem, {0: 0, 1: 1}), RefinementConfig(event_slack=0))
    candidate = decoded.candidates[0]
    assert query_assignments(decoded, candidate, {0: 0}).status == 'allowed'

    def forbidden(*args, **kwargs):
        raise AssertionError('An impossible mapping needs no event query')
    monkeypatch.setattr('mappa.mappa_core.core.symmetry.event_counts._event_model', forbidden)
    assert query_assignments(decoded, candidate, {0: 2}).status == 'forbidden'


def test_empty_mapping_has_no_events_and_creates_no_atoms():
    problem=AAMProblem(endpoint('CO',[(0,1,1)]),endpoint('CCC',[(0,1,1),(1,2,1)]))
    decoded=refine_mapping(catalogue(problem,{}),RefinementConfig(event_slack=0))
    assert decoded.complete and decoded.minimum_proven and decoded.observed_minimum==0
    assert decoded.candidates[0].mapping=={}


@pytest.mark.parametrize('joined', [False, True])
def test_correlated_local_choices_preserve_partial_event_union(joined):
    problem = AAMProblem(endpoint('CCCC', [(0,1,1), (1,2,1), (2,3,1)]),
                         endpoint('CCCCC', [(0,1,1), (1,2,.6), (2,3,1.8), (0,3,.4)]))
    first = saved_fragment(problem, source_atoms=(0,1), target_atoms=(0,1))
    continuation_atoms = (0,1,2) if joined else (1,2)
    anchors = ((0,0), (1,1)) if joined else ((1,1),)
    last = saved_fragment(problem, source_atoms=continuation_atoms,
                          target_atoms=continuation_atoms, anchors=anchors)
    mapping = {0:0, 1:1, 2:2}
    family = FinalFamily(tuple(mapping.items()), (),
                         (('relation',last), ('relation',first)), (.2,1.))
    # Independently compose each complete local relation in recorded order.
    mappings = {tuple(mapping.items())}
    for fragment in (last, first):
        expanded = set()
        for pairs in represented((fragment,)):
            chosen = dict(pairs)
            transport = {p:chosen[r] for r,p in fragment.witness}
            expanded.update(tuple((r, transport.get(p,p)) for r,p in old) for old in mappings)
        mappings = expanded
    assert len(mappings) == 2
    index = SignedEventIndex(problem)
    expected = {index.describe(index.vector(dict(pairs)))['id'] for pairs in mappings}
    decoded = extract_path_events(family.as_path(problem), problem, index,
                                  seconds=30, max_patterns=None)
    assert decoded['complete']
    assert {pattern['id'] for pattern in decoded['patterns']} == expected


def test_distinct_fixed_boundary_elements_match_exhaustive_events():
    # N and O each have one possible image: both are plain integer constants.
    problem = AAMProblem(endpoint('C', []),
                         endpoint('CCNO', [(0, 2, 1), (0, 3, 1), (1, 2, 1)]))
    cat = catalogue(problem, {0: 0}, (0, 1))
    index = SignedEventIndex(problem)
    expected = {index.describe(index.vector({0: p}))['id'] for p in (0, 1)}
    result = refine_mapping(cat)
    assert result.complete and result.minimum_proven and result.observed_minimum == 1
    assert {candidate.id for candidate in result.candidates} == expected
    assert sorted(candidate.total for candidate in result.candidates) == [1, 2]
    for candidate in result.candidates:
        assert candidate.total == sum(oracle(problem, candidate.mapping))
