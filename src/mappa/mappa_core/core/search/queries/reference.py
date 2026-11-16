"""Shared fragment membership before composing an exact reference witness."""
from collections import OrderedDict
import time

from mappa.mappa_core.core.symmetry.actions import ActionIndex, _PreparedAction
from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
from .family import compile_actions, realize_actions


class ReferenceQuery:
    """One reference, sharing local answers across saved branch compositions.

    Disjoint action supports on annotated element classes are independent.
    Overlapping supports keep their complete ordered joint query. Unannotated
    atoms remain in every local relation, so a heavy-atom answer still requires
    a complete injective discovery. Every action is an element-preserving
    permutation: composing its full witness also preserves unannotated occupancy.
    """
    def __init__(self, problem, reference, *, timeout_ms=1000):
        self.problem, self.reference = problem, dict(reference)
        self.timeout_ms = timeout_ms
        self.index = ActionIndex(problem.target_atom_count)
        self.elements = problem.product.elements
        self.annotated_elements = {problem.reactant.elements[r] for r in reference}
        self.actions, self.supports = {}, {}
        self.components, self.answers = OrderedDict(), OrderedDict()
        self.stats = dict(fragment_queries=0, joint_queries=0, cache_hits=0,
                          composition_rejections=0, direct_witnesses=0)
        self.valid = (len(set(reference.values())) == len(reference) and all(
            0 <= p < len(self.elements) and problem.reactant.elements[r] == self.elements[p]
            for r, p in reference.items()))

    def _prepare(self, action):
        action = _PreparedAction(action) if not isinstance(action, _PreparedAction) else action
        if action not in self.actions:
            bounds = self.index.action_domains(action)
            # This is a property of the saved relation, not of one chosen witness.
            if any(self.elements[p] != self.elements[q]
                   for p, domain in enumerate(bounds) for q in domain):
                raise ValueError('Reference composition requires element-preserving saved actions')
            self.actions[action] = action
            self.supports[action] = frozenset(p for p, domain in enumerate(bounds)
                if domain != frozenset((p,)) and self.elements[p] in self.annotated_elements)
        return self.actions[action]

    @staticmethod
    def _remember(cache, key, value):
        cache[key] = value
        cache.move_to_end(key)
        if len(cache) > 4096:
            cache.popitem(last=False)
        return value

    def _components(self, actions):
        if actions not in self.components:
            groups = []
            for i, action in enumerate(actions):
                support = self.supports[action]
                if not support:
                    continue
                positions = [i]
                separate = []
                for old_support, old_positions in groups:
                    if support & old_support:
                        support = support | old_support
                        positions += old_positions
                    else:
                        separate.append((old_support, old_positions))
                separate.append((support, positions))
                groups = separate
            self._remember(self.components, actions, tuple(
                (support, tuple(sorted(positions))) for support, positions in groups))
        return self.components[actions]

    def query(self, mapping, program):
        """Test a composition; never expand its individual bijections."""
        mapping = dict(mapping)
        if not self.valid or any(r not in mapping for r in self.reference):
            return 'not_recovered', None
        actions = tuple(self._prepare(action) for action in program)
        conditions = {mapping[r]: p for r, p in self.reference.items()}
        pending, covered = [], set()
        for support, positions in self._components(actions):
            required = tuple(sorted((p, q) for p, q in conditions.items() if p in support))
            covered.update(p for p, _ in required)
            if not required:
                continue
            if any(q not in support for _, q in required):
                self.stats['composition_rejections'] += 1
                return 'not_recovered', None
            key = tuple(actions[i] for i in positions), required
            if key in self.answers and self.answers[key][0] == 'not_recovered':
                self.stats['cache_hits'] += 1
                return 'not_recovered', None
            pending.append((key, positions))
        if any(p != q for p, q in conditions.items() if p not in covered):
            self.stats['composition_rejections'] += 1
            return 'not_recovered', None
        selected, unknown = {}, False
        # Check local fragments before any genuinely coupled placement program.
        for key, positions in sorted(pending, key=lambda item: len(item[0][0])):
            if key in self.answers:
                self.stats['cache_hits'] += 1
                status, choices = self.answers[key]
            else:
                status, choices = self._query_component(*key)
                self._remember(self.answers, key, (status, choices))
            if status == 'not_recovered':
                return status, None
            unknown |= status == 'unknown'
            if status == 'recovered':
                selected.update((positions[i], choice) for i, choice in choices.items())
        if unknown:
            return 'unknown', None
        actual = mapping
        evidence = []
        for i, action in enumerate(actions):
            choice = selected.get(i)
            if choice is None:
                continue  # The saved witness is the identity action.
            transport = dict(choice['permutation'])
            actual = {r: transport.get(p, p) for r, p in actual.items()}
            evidence.append(dict(choice, action=i))
        if (len(set(actual.values())) != len(actual)
                or any(actual[r] != p for r, p in self.reference.items())
                or any(self.problem.reactant.elements[r] != self.elements[p]
                       for r, p in actual.items())):
            raise ValueError('Composed reference certificate is inconsistent')
        return 'recovered', dict(mapping=sorted(actual.items()), actions=evidence,
            scope='full_explicit', certificate='shared_fragment_composition')

    def _query_component(self, actions, required):
        if all(p == q for p, q in required):
            return 'recovered', {}
        if len(actions) == 1 and actions[0][0] == 'relation':
            self.stats['fragment_queries'] += 1
            direct = self._direct_fragment(actions[0][1], required)
            if direct is not None:
                self.stats['direct_witnesses'] += 1
                return 'recovered', {0: direct}
        else:
            self.stats['joint_queries'] += 1
        # Atom domains reject impossible requests but never prove acceptance.
        for p, wanted in required:
            domain = frozenset((p,))
            for action in actions:
                bounds = self.index.action_domains(action)
                domain = frozenset(q for p in domain for q in bounds[p])
            if wanted not in domain:
                return 'not_recovered', None
        import z3
        started = time.perf_counter()
        encoder = FiniteConstraints()
        values, compiled = compile_actions(encoder,
            [encoder.constant(p) for p, _ in required],
            [(kind, i, data) for i, (kind, data) in enumerate(actions)])
        for value, (_, wanted) in zip(values, required):
            encoder.restrict(value[0], {wanted})
        encoder.propagate()
        if encoder.failed:
            return 'not_recovered', None
        left = self.timeout_ms - 1000 * (time.perf_counter() - started)
        if left <= 0:
            return 'unknown', None
        status, model, _ = encoder.solve(max(1, int(left)))
        if status != z3.sat:
            return ('not_recovered' if status == z3.unsat else 'unknown'), None
        choices = realize_actions(encoder, compiled, model, len(self.elements))
        return 'recovered', {i: choice for i, choice in enumerate(choices)}

    def _direct_fragment(self, fragment, required):
        """Validate a supplied full discovery before constructing constraints."""
        relation = fragment.relation
        inverse = {p: r for r, p in fragment.witness}
        discovery = dict(relation.anchors)
        for p, q in required:
            r = inverse[p]
            if r in discovery and discovery[r] != q:
                return None
            discovery[r] = q
        if len(set(discovery.values())) != len(discovery):
            return None
        unused = set(fragment.target_atoms) - set(discovery.values())
        for r, p in fragment.witness:
            if r not in discovery and p in unused:
                discovery[r] = p
                unused.remove(p)
        for r, domain in relation.domains:
            if r not in discovery:
                p = next((p for p in domain if p in unused), None)
                if p is None:
                    return None
                discovery[r] = p
                unused.remove(p)
        if not relation.allows_discovery(discovery):
            return None
        identity = tuple(range(len(fragment.source_atoms)))
        return dict(kind='relation', transition=0, discovery=discovery,
            source_permutation=identity, target_permutation=identity,
            permutation=sorted((p, discovery[r]) for r, p in fragment.witness))
