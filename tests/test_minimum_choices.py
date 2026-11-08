"""Independent exhaustive checks of the minimum-choice experiment on small graphs."""
from dataclasses import replace
from pathlib import Path
import sys
import time

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from minimum_choices import ProjectedChoices, chemical_orbit, count_choices, heavy_images
from slap_audit import equivalent, paired_graph
from mappa import AAMProblem, decode_events
from growth_cases import saved_fragment
from test_event_interpretation import endpoint, saved


def feature(ep):
    heavy = [i for i,e in enumerate(ep.elements) if e != 'H']
    return dict(heavy=heavy, colors=[(ep.elements[i],) for i in heavy],
        bonds=[(i,j,(float(ep.wbo[a,b]),)) for i,a in enumerate(heavy)
               for j,b in enumerate(heavy) if i < j and ep.wbo[a,b] > .2])


def all_operations(group):
    """Independent closure; tiny test graphs only."""
    identity = tuple(range(group.degree))
    found, pending = {identity}, [identity]
    while pending:
        old = pending.pop()
        for g in group.generators:
            new = tuple(old[g[i]] for i in range(group.degree))
            if new not in found:
                found.add(new); pending.append(new)
    return found


@pytest.mark.parametrize('source_edges,target_edges', [
    ([(0,1,1.),(1,2,2.)], [(0,1,2.),(1,2,1.)]),
    ([(0,1,1.),(1,2,1.),(2,0,1.)], [(0,1,1.),(1,2,1.)]),
    ([(0,1,1.),(1,2,1.)], [(0,1,1.),(1,2,1.),(2,0,1.)]),
    ([(0,1,1.),(1,2,2.),(2,3,1.)], [(0,1,1.),(0,2,2.),(0,3,1.)]),
])
def test_every_saved_event_count_matches_exhaustive_joint_actions(source_edges, target_edges):
    n = 1+max(i for edges in (source_edges,target_edges) for a,b,_ in edges for i in (a,b))
    r,p = endpoint('C'*n, source_edges), endpoint('C'*n, target_edges)
    fragment = saved_fragment(AAMProblem(r,p), iso_tolerance=3.)
    decoded = decode_events(saved(r,p,dict(enumerate(range(n))),fragments=(fragment,)))
    table = decoded.tables(0)
    oracle = {}
    for operation in all_operations(table.program.group):
        mapping = table.program.mapping(np.asarray(operation))
        score = table.evaluate(decoded.index.vector(mapping))['total']
        oracle.setdefault(score,set()).add(tuple(sorted(mapping.items())))
    assert heavy_images(table.program,set(range(n)),deadline=time.perf_counter()+5) == set().union(*oracle.values())
    for score in range(max(oracle)+2):
        projected = ProjectedChoices(table, score, time.perf_counter()+10.)
        found = set()
        while (value := projected.next(set(range(n)))) is not None:
            assert value not in found
            found.add(value)
        assert found == oracle.get(score,set())


def test_count_is_mapping_equivalence_not_branches_directions_or_hydrogens():
    r = endpoint('CCHHHH', [(0,1,1.),(0,2,1.),(0,3,1.),(1,4,1.),(1,5,1.)])
    f = saved_fragment(AAMProblem(r,r))
    catalogue = saved(r,r,dict(enumerate(range(6))),fragments=(f,))
    catalogue.families.append(catalogue.families[0])
    next(iter(catalogue.branches.values())).add(1)
    decoded = decode_events(catalogue)
    reference = dict(features=[feature(r)]*2,mapping=[(0,0),(1,1)])
    result = count_choices([(False,decoded),(True,decoded)],reference)
    assert result['status'] == 'complete'
    assert result['equivalent_labels_excluded'] == 2
    assert result['choices'] == 1 and result['random_success'] == 1
    table = decoded.tables(0)
    assert heavy_images(table.program,{0,1},deadline=time.perf_counter()+5) == {
        ((0,0),(1,1)),((0,1),(1,0))}


def test_reference_in_branch_is_not_reference_at_minimum():
    r = endpoint('CCC',[(0,1,1.),(1,2,2.)])
    p = endpoint('CCC',[(0,1,2.),(1,2,1.)])
    f = saved_fragment(AAMProblem(r,p),iso_tolerance=2.)
    decoded = decode_events(saved(r,p,{0:0,1:1,2:2},fragments=(f,)))
    reference = dict(features=[feature(r),feature(p)],mapping=[(0,0),(1,1),(2,2)])
    result = count_choices([(False,decoded)],reference)
    assert result['choices'] == 1 and result['detected_minimum'] == 0
    assert result['reference_at_minimum'] is False and result['random_success'] == 0


def test_exhausted_count_is_not_an_exact_small_shortlist():
    r = endpoint('CC',[(0,1,1.)])
    f = saved_fragment(AAMProblem(r,r))
    decoded = decode_events(saved(r,r,{0:0,1:1},fragments=(f,)))
    reference = dict(features=[feature(r)]*2,mapping=[(0,0),(1,1)])
    result = count_choices([(False,decoded)],reference,max_mappings=1)
    assert result['status'] == 'incomplete' and result['choices'] is None
    assert result['choices_lower_bound'] == 1 and result['random_success'] is None


def test_multiple_real_options_and_partial_mapping_equivalence():
    r = endpoint('CCCO',[(0,1,1.),(1,2,1.),(2,3,1.)])
    p = endpoint('CCCO',[(0,1,1.),(1,2,1.),(2,3,1.)])
    # Fixed mappings have no local search; different family partitions keep their
    # existing model count. Independent VF2 is the class-count oracle.
    catalogue = saved(r,p,{0:0,1:1})
    for mapping in ({0:1,1:2},{0:2,1:1},{0:1,1:0}):
        other = saved(r,p,mapping)
        number = len(catalogue.families)
        catalogue.families.extend(other.families)
        for key in other.branches:
            catalogue.branches.setdefault(key,set()).add(number)
    decoded = decode_events(catalogue)
    reference = dict(features=[feature(r),feature(p)],mapping=[(0,0),(1,1)])
    expected = []
    for context in decoded.window.contexts:
        table = decoded.tables(*context[:2],cut=context[2])
        pairs = table.program.base
        if table.evaluate(decoded.index.vector(dict(pairs)))['total'] == decoded.minimum_events:
            graph = paired_graph(reference['features'],pairs)
            if not any(equivalent(graph,old) is not None for old in expected):
                expected.append(graph)
    result = count_choices([(False,decoded)],reference)
    assert result['choices'] == len(expected) > 1
    assert result['random_success'] == 1/len(expected)


@pytest.mark.parametrize('placement', [False, True])
def test_sweep_and_placement_counts_preserve_joint_actions(placement):
    if placement:
        from test_final_branches import placed_fragment_path
        from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
        problem,path = placed_fragment_path()
        catalogue = FinalBranchCatalogue(problem).add_graph(path.graph)
    else:
        r = endpoint('CCCC',[(0,1,1.),(1,3,1.)])
        p = endpoint('CCCC',[(0,2,1.),(1,3,1.)])
        f = saved_fragment(AAMProblem(r,p),source_atoms=(1,2,3),anchors=((3,3),))
        catalogue = saved(r,p,dict(enumerate(range(4))),fragments=(f,),
                          pairs=(((0,),(0,)),((1,2,3),(1,2,3))))
        catalogue.families[0] = replace(catalogue.families[0],provenance=[dict(cuts=((0,1),))])
    decoded = decode_events(catalogue)
    for context in decoded.window.contexts:
        table = decoded.tables(*context[:2],cut=context[2])
        expected = set()
        for operation in all_operations(table.program.group):
            mapping = table.program.mapping(np.asarray(operation))
            if table.evaluate(decoded.index.vector(mapping))['total'] == decoded.minimum_events:
                expected.add(tuple(sorted(mapping.items())))
        all_images = {tuple(sorted(table.program.mapping(np.asarray(operation)).items()))
                      for operation in all_operations(table.program.group)}
        assert heavy_images(table.program,set(dict(table.program.base)),deadline=time.perf_counter()+5) == all_images
        projected = ProjectedChoices(table,decoded.minimum_events,time.perf_counter()+10)
        observed = set()
        while (pairs := projected.next(set(dict(table.program.base)))) is not None:
            observed.add(pairs)
        assert observed == expected


def test_orbit_exclusions_are_only_chemical_equivalence_including_unmatched_atoms():
    from golden import endpoint_generators
    r = endpoint('CCCC',[(0,1,1.),(1,2,1.),(2,3,1.)])
    features = [feature(r)]*2
    generators = [endpoint_generators(f,4) for f in features]
    for pairs in (((0,0),(1,1)),((0,0),(1,2),(2,1),(3,3))):
        graph = paired_graph(features,pairs)
        images = chemical_orbit(pairs,generators)
        assert len(images) > 1
        assert all(equivalent(graph,paired_graph(features,image)) is not None for image in images)
        assert chemical_orbit(pairs,generators,limit=1) == {pairs}
