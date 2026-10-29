from mappa.mappa_workflow.mapping_analytics.queries import query_assignments
"""One event rule for partial and full maps; exhaustive tiny graphs are the oracle."""
from itertools import combinations, permutations
import random

import numpy as np
import pytest

from mappa import AAMProblem, MolecularEndpoint
from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_workflow.mapping_analytics.event_query import _event_model, _same_event_orbit
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue, FinalFamily
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from mappa.mappa_core.post_processing.decoder.api import decode_events
from mappa.mappa_workflow.mapping_analytics.family import compile_path
from mappa.mappa_core.post_processing.decoder.events import classify_bonds
from reference_event_decoder import extract_path_events


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
            count += source.wbo[a, b] - other >= threshold
        result.append(int(count))
    return tuple(result)


def catalogue(problem, mapping, pool=()):
    cat = FinalBranchCatalogue(problem)
    cat.families = [FinalFamily(tuple(sorted(mapping.items())), (),
        (('pool', tuple(pool)),) if pool else (), (.2, 1.))]
    pairs = ((tuple(sorted(mapping)), tuple(sorted(mapping.values()))),) if mapping else ()
    cat.branches[((.2, 1.), pairs)] = {0}
    return cat


def test_disconnected_unmatched_material_is_ignored_in_both_directions():
    core = endpoint('C', [])
    with_water = endpoint('COHH', [(1, 2, 1), (1, 3, 1)])
    for r, p in [(core, with_water), (with_water, core)]:
        problem = AAMProblem(r, p)
        decoded = decode_events(catalogue(problem, {0: 0}), EventDecodeConfig(event_slack=0))
        assert decoded.complete and decoded.minimum_proven and decoded.observed_minimum == 0
        assert decoded.candidates[0].mapping == {0: 0}
        assert decoded.candidates[0].endpoint_events == {'broken': (), 'formed': ()}


def test_mapped_boundary_and_interfragment_bonds_use_same_rule():
    r = endpoint('CCOH', [(0,1,1), (1,2,1), (2,3,1)])
    p = endpoint('CCNH', [(0,2,1), (2,3,1)])
    problem = AAMProblem(r, p)
    decoded = decode_events(catalogue(problem, {0:0,1:1}), EventDecodeConfig(event_slack=0))
    assert decoded.minimum_proven and decoded.observed_minimum == 3
    candidate = decoded.candidates[0]
    assert candidate.endpoint_events == {'broken': ((0,1),(1,2)), 'formed': ((0,2),)}
    reverse = decode_events(catalogue(AAMProblem(p,r), {0:0,1:1}), EventDecodeConfig(event_slack=0))
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
        path = cat.families[0].as_path(problem)
        compiled = compile_path(path,problem,{},source_atoms=(),complete_reference=False)
        terms,total,_ = _event_model(compiled,index,time.perf_counter()+30)
        equalities = [_same_event_orbit(terms,total,p,index) for p in patterns[:3]]
        for mapping, pattern in zip(maps,patterns):
            compiled.solver.push()
            compiled.solver.add(*(compiled.values[a][0]==p for a,p in mapping.items()))
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
        assert {p['id'] for p in result['patterns']} == {p['id'] for p in patterns}


def test_partial_target_symmetry_needs_no_solver(monkeypatch):
    r = endpoint('C', [])
    p = endpoint('CCC', [(0,1,1),(1,2,1),(0,2,1)])
    cat = catalogue(AAMProblem(r,p),{0:0},range(3))
    def forbidden(*args,**kwargs):
        raise AssertionError('Exact target symmetry needs no solver')
    monkeypatch.setattr('mappa.mappa_workflow.mapping_analytics.family.compile_path',forbidden)
    result=decode_events(cat,EventDecodeConfig(event_slack=0))
    assert result.minimum_proven and result.observed_minimum==2
    assert len(result.candidates)==1 and result.family_reports[0]['solver_queries']==0


def test_saved_partial_graph_survives_catalogue_and_artifact():
    from mappa.mappa_core.utils.io.records import decoded_events_record
    problem=AAMProblem(endpoint('CO',[(0,1,1)]),endpoint('CON',[(0,2,1)]))
    family=catalogue(problem,{0:0,1:1}).families[0]
    cat=FinalBranchCatalogue(problem).add_graph(family.as_path(problem).graph)
    assert len(cat.families)==1 and cat.families[0].mapping==family.mapping
    decoded=decode_events(cat,EventDecodeConfig(event_slack=0))
    record=decoded_events_record(decoded)
    assert record['complete'] and record['observed_minimum']==2
    assert record['candidates'][0]['endpoint_events']['formed']==((0,2),)
    candidate=decoded.candidates[0]
    assert query_assignments(decoded, candidate,{0:0,1:1}).status=='allowed'


def test_empty_mapping_has_no_events_and_creates_no_atoms():
    problem=AAMProblem(endpoint('CO',[(0,1,1)]),endpoint('CCC',[(0,1,1),(1,2,1)]))
    decoded=decode_events(catalogue(problem,{}),EventDecodeConfig(event_slack=0))
    assert decoded.complete and decoded.minimum_proven and decoded.observed_minimum==0
    assert decoded.candidates[0].mapping=={}


@pytest.mark.parametrize('actions', [
    (('group', ((1,2,3,0),)),),
    (('pool', (0,1,2)), ('pool', (1,2,3))),
    (('group', ((1,2,3,0),)), ('pool', (0,1))),
])
def test_correlated_and_overlapping_actions_preserve_partial_event_union(actions):
    problem=AAMProblem(endpoint('CCC',[(0,1,1),(1,2,1)]),
                       endpoint('CCCC',[(0,1,1),(1,2,.6),(2,3,1.8),(0,3,.4)]))
    mapping={0:0,2:1}
    cat=catalogue(problem,mapping)
    cat.families[0].actions=actions
    vectors={tuple(mapping[a] for a in sorted(mapping))}
    for kind,data in actions:
        if kind=='pool':
            permutations_=[]
            for order in permutations(data):
                g=list(range(4))
                for a,b in zip(data,order):g[a]=b
                permutations_.append(g)
        else:
            group={tuple(range(4))}
            pending=list(group)
            while pending:
                old=pending.pop()
                for generator in data:
                    new=tuple(generator[p] for p in old)
                    if new not in group:group.add(new);pending.append(new)
            permutations_=group
        vectors={tuple(g[p] for p in v) for v in vectors for g in permutations_}
    index=SignedEventIndex(problem)
    expected={index.describe(index.vector(dict(zip(sorted(mapping),v))))['id'] for v in vectors}
    decoded=extract_path_events(cat.families[0].as_path(problem), problem, index, seconds=30, max_patterns=None)
    assert decoded['complete']
    assert {c['id'] for c in decoded['patterns']}==expected


def test_distinct_fixed_boundary_elements_match_exhaustive_events():
    # N and O each have one possible image: both are plain integer constants.
    problem = AAMProblem(endpoint('C', []),
                         endpoint('CCNO', [(0, 2, 1), (0, 3, 1), (1, 2, 1)]))
    cat = catalogue(problem, {0: 0}, (0, 1))
    index = SignedEventIndex(problem)
    expected = {index.describe(index.vector({0: p}))['id'] for p in (0, 1)}
    result = decode_events(cat)
    # Conversion keeps this saved singleton target set; it does not enumerate
    # the old pool's alternative target regions.
    assert result.complete and result.minimum_proven and result.observed_minimum == 2
    assert {candidate.id for candidate in result.candidates} <= expected
    assert [candidate.total for candidate in result.candidates] == [2]
    for candidate in result.candidates:
        assert candidate.total == sum(oracle(problem, candidate.mapping))
