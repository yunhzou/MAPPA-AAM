"""Exact symbolic queries of saved decoded choices, shared by benchmark analyses.

This external analysis encodes the saved ordered operations and event model.
It does not search for new fragment starts or alter the decoded response.
"""
from types import SimpleNamespace
import time

import numpy as np
import z3

from mappa.mappa_core.core.symmetry.constraints import SymbolicActions


class ProjectedChoices:
    """Project the saved ordered operations onto mappings and their saved score."""
    def __init__(self, table, event_count, deadline, *, heavy=None, known=(), assignments=None):
        self.table, self.deadline = table, deadline
        self.solver = z3.Solver()
        self.encoder = SymbolicActions(self.solver, fresh=True)
        program = table.program
        self.prefix = self.encoder.act([self.encoder.constant(i) for i in range(program.size)],
                                      program.group.generators, 'saved_decoded_choices')
        self.values = self.mapping(program, self.prefix)
        self.assignments = dict(assignments or {})
        for atom, image in self.assignments.items():
            self.solver.add(self.values[atom][0] == image if atom in self.values else False)
        self.domain_solver = z3.Solver()
        self.domain_values = {a: z3.FreshInt('heavy_image') for a in self.values
                              if heavy is None or a in heavy}
        for a,value in self.domain_values.items():
            self.domain_solver.add(z3.Or(*(value == p for p in self.values[a][1])))
        self.domain_solver.add(z3.Distinct(*self.domain_values.values()))
        for pairs in known:
            if set(dict(pairs)) == set(self.domain_values) and all(b in self.values[a][1] for a,b in pairs):
                self.exclude(pairs)
        self.event_count = z3.FreshInt('event_count') if isinstance(event_count, tuple) else event_count
        if isinstance(event_count, tuple):
            self.solver.add(self.event_count >= event_count[0], self.event_count <= event_count[1])
        self.covered = self.domain_covered()
        if self.covered:
            return
        constant = sum(int(np.asarray(row).sum()) for row in table.fixed_boundary.values())
        terms = []
        for group in table._term_groups():
            for term in group:
                self.check()
                response = term.response
                if response.constant is not None:
                    constant += int(bool(response.constant))
                    continue
                a, b = (self.event_value(i) for i in term.atoms)
                condition = self.pair_condition(a, b,
                    lambda x, y: bool(response.lookup(table.index.p[x, y])))
                terms.append(z3.If(condition, 1, 0))
        if table.sweep_cut is not None:
            constant -= int(table.boundary_events.source_cut_counted)
            a, b = (self.event_value(i) for i in table.sweep_cut)
            condition = self.pair_condition(a, b,
                lambda x, y: tuple(sorted((x, y))) in table.product_boundary_set)
            terms.append(-z3.If(condition, 1, 0))
        self.solver.add(constant + z3.Sum(terms) == self.event_count)

    def event_value(self, atom):
        # The exact requested assignments make most heavy-atom costs constant.
        return (self.encoder.constant(self.assignments[atom]) if atom in self.assignments
                else self.values[atom])

    def exclude(self, pairs):
        self.solver.add(z3.Or(*(self.values[a][0] != b for a,b in pairs)))
        self.domain_solver.add(z3.Or(*(self.domain_values[a] != b for a,b in pairs)))

    def domain_covered(self):
        """A relaxed injective image domain can certify no unseen heavy output.

        This ignores event and symmetry restrictions, so only UNSAT can skip the
        actual query. It avoids proving absence through irrelevant H operations.
        """
        self.check()
        self.domain_solver.set(timeout=max(1,min(50,int(1000*(self.deadline-time.perf_counter())))))
        return self.domain_solver.check() == z3.unsat

    def check(self):
        if time.perf_counter() >= self.deadline:
            raise TimeoutError('choice-count allowance exhausted')

    def remap(self, value, transform):
        return self.encoder.lookup(value, {i: self.encoder.constant(transform(i)) for i in value[1]})

    def mapping(self, response, prefix):
        self.check()
        if hasattr(response, 'actions'):
            values = {a: self.encoder.constant(b) for a, b in response.base}
            for ((kind, data), prepared), (start, atoms, degree) in zip(response.actions, response.frames):
                local = [self.remap(value, lambda p: p-start) for value in prefix[start:start+degree]]
                if kind == 'relation':
                    selected = self.mapping(prepared, local)
                    transport = {b: selected[a] for a, b in data.witness}
                else:
                    transport = {a: self.remap(value, lambda p: atoms[p])
                                 for a, value in zip(atoms, local)}
                values = {a: self.encoder.lookup(value, transport) for a, value in values.items()}
            return values
        n = len(response.witness)
        targets = [self.remap(value, lambda p: response.fragment.target_atoms[response.witness[p]])
                   for value in (prefix if not response.two_sided else
                                 [self.remap(v, lambda p: p-n) for v in prefix[n:]])]
        if response.two_sided:
            # Source and target operations are joint permutations, never free pools.
            inverse = []
            for atom in range(n):
                possible = [i for i, value in enumerate(prefix[:n]) if atom in value[1]]
                if len(possible) == 1:
                    inverse.append(self.encoder.constant(possible[0]))
                else:
                    value = self.encoder.symbol(possible)
                    self.solver.add(z3.Or(*(z3.And(value[0] == i, prefix[i][0] == atom) for i in possible)))
                    inverse.append(value)
            targets = [self.encoder.lookup(value, dict(enumerate(targets))) for value in inverse]
        return dict(zip(response.fragment.source_atoms, targets))

    def pair_condition(self, left, right, predicate):
        rows = []
        for a in left[1]:
            self.check()
            targets = [b for b in right[1] if predicate(a, b)]
            if targets:
                rows.append(z3.And(left[0] == a, z3.Or(*(right[0] == b for b in targets))))
        return z3.simplify(z3.Or(*rows))

    def find(self, *, max_conflicts=None):
        self.check()
        if self.covered or self.domain_covered():
            return None
        self.solver.set(timeout=max(1, int(1000*(self.deadline-time.perf_counter()))))
        if max_conflicts is not None:
            self.solver.set(max_conflicts=max_conflicts)
        status = self.solver.check()
        if status == z3.unknown:
            raise TimeoutError(self.solver.reason_unknown())
        if status == z3.unsat:
            return None
        model = self.solver.model()
        def number(value):
            return value if isinstance(value, int) else model.eval(value, model_completion=True).as_long()
        prefix = np.asarray([number(value[0]) for value in self.prefix])
        count = number(self.event_count)
        mapping = self.table.program.mapping(prefix)
        if mapping != {a: number(value[0]) for a, value in self.values.items()}:
            raise AssertionError('Symbolic mapping differs from saved ordered action replay')
        vector = self.table.index.vector(mapping)
        if (self.table.program.fixed_events + self.table.program.score(prefix, self.table.program.pairs)
                != count or self.table.evaluate(vector)['total'] != count):
            raise AssertionError('Choice does not have the saved requested event count')
        self.table.index._vectors([vector])
        witness, choices = self.table.program.query_result(prefix, count, 0., None)
        return SimpleNamespace(status='found', witness=witness, choices=choices, event_count=count)

    def next(self, heavy):
        answer = self.find()
        if answer is None:
            return None
        mapping = dict(answer.witness)
        pairs = tuple(sorted((a, b) for a, b in mapping.items() if a in heavy))
        self.exclude(pairs)
        return pairs

