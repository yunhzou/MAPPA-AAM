"""Exactness checks for performance reductions, without permutation enumeration."""
from itertools import combinations
import time
from types import SimpleNamespace

import numpy as np
import z3

from mappa.mappa_core.post_processing.decoder.events import SignedEventIndex
from mappa.mappa_workflow.mapping_analytics.event_query import _event_model, _event_count_constraint
from mappa.mappa_workflow.mapping_analytics.family import compile_path
from mappa.mappa_core.post_processing.final_branches import FinalFamily
from mappa.mappa_core.post_processing.decoder.events import classify_bonds
from mappa.mappa_core.utils.bonds import is_metal_element
from reference_event_model import _event_model as reference_event_model
from test_event_patterns import problem


def test_permutation_composition_implies_injectivity_with_overlapping_pools():
    value = problem('OOOOO', [], [])
    family = FinalFamily(tuple(enumerate(range(5))), (),
                         (('pool', (0, 1, 2)), ('group', ((1, 0, 3, 2, 4),)),
                          ('pool', (2, 3, 4))), (.5, 1.))
    compiled = compile_path(family.as_path(value), value, {}, source_atoms=(), complete_reference=False)
    images = [v[0] for v in compiled.values]
    # No global Distinct premise: prove it from the correlated local actions.
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
        # Deliberately allow broad supports; element restrictions need not be
        # assumed by the local table cache to preserve exact signed responses.
        values = [(z3.Int(f'image_{i}'), frozenset(range(4))) for i in range(4)]
        compiled = SimpleNamespace(values=values, solver=z3.Solver(), representative=dict(enumerate(range(4))))
        compiled.solver.add(z3.Distinct(*(v[0] for v in values)),
                            *(z3.And(v[0] >= 0, v[0] < 4) for v in values))
        old, total, metrics = reference_event_model(compiled, index, time.perf_counter()+5)
        for _ in range(2):  # Cold build followed by cache reuse.
            new, new_total, new_metrics = _event_model(compiled, index, time.perf_counter()+5)
            assert new_metrics == metrics
            differences = [total != new_total] + [old[k] != new[k] for k in old]
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
    from mappa.mappa_core.post_processing.decoder.events import FragmentEventCounter
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
    from mappa.mappa_core.post_processing.decoder.api import decode_events
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
    result = decode_events(source)
    assert result.candidates[0].id == expected['id']
    assert result.candidates[0].endpoint_events == expected['endpoint_events']
    assert result.candidates[0].total == expected['total']
    assert result.counting_statistics['branch_counts'] == len(result.family_reports)


def test_many_combinations_reuse_fragment_and_boundary_counts():
    from itertools import product
    from mappa.mappa_core.post_processing.decoder.events import FragmentEventCounter
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
