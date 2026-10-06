"""Exactness checks for performance reductions, without permutation enumeration."""
from itertools import combinations
import time
from types import SimpleNamespace

import numpy as np
import z3

from mappa.event_patterns import SignedEventIndex, _event_model, _event_count_constraint
from mappa.family_query import compile_path
from mappa.final_branches import FinalFamily
from mappa.frag import classify_bonds, is_metal_element
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
        compiled = SimpleNamespace(values=values, solver=z3.Solver())
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


def test_sparse_scalar_validation_matches_original_dense_rule():
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
                    tau = metal if metal is not None and any(is_metal_element(elements[i]) for i in (a,b)) else .5
                    w=source[a,b]
                    if w < tau:continue
                    mapped = a in m and b in m
                    other=target[m[a],m[b]] if mapped else None
                    if not mapped or w-other >= tau:out.append((a,b,float(w),None if other is None else float(other)))
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
