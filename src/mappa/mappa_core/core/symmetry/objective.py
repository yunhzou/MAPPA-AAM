"""Original-weight event counts and a connectivity-free lower bound."""
from collections import defaultdict
from itertools import combinations

import numpy as np

from mappa.mappa_core.utils.bonds import is_metal_element


class FragmentObjective:
    def __init__(self, fragment, event_problem=None, threshold=.5, metal_threshold=.3):
        source, target = fragment.source_symmetry, fragment.target_symmetry
        if event_problem is None:
            self.source, self.target = source.weights, target.weights
            self.source_elements, self.target_elements = source.elements, target.elements
        else:
            left, right = event_problem.reactant, event_problem.product
            self.source = self._weights(left, source)
            self.target = self._weights(right, target)
            self.source_elements = tuple(left.elements[a] for a in source.atoms)
            self.target_elements = tuple(right.elements[a] for a in target.atoms)
        self.size = len(source.atoms)
        self.thresholds = np.broadcast_to(float(threshold), (self.size, self.size))
        if metal_threshold is not None:
            metal = np.asarray([is_metal_element(element) for element in self.source_elements], dtype=bool)
            if np.any(metal) and metal_threshold != threshold:
                self.thresholds = self.thresholds.copy()
                self.thresholds[metal[:, None] | metal[None, :]] = float(metal_threshold)
        self.thresholds.setflags(write=False)
        self.target_positions = {atom: i for i, atom in enumerate(target.atoms)}
        self.source_atoms = source.atoms
        self.anchors = tuple((source.atoms.index(a), target.atoms.index(b))
                             for a, b in fragment.fixed_assignments)

    @staticmethod
    def _weights(endpoint, symmetry):
        weights = endpoint.wbo[np.ix_(symmetry.atoms, symmetry.atoms)]
        # Discovery and event weights may differ. Share only an exact match;
        # keep the caller's original event weights otherwise.
        if np.array_equal(weights, symmetry.weights):
            return symmetry.weights
        weights.setflags(write=False)
        return weights

    def images(self, witness):
        mapping = dict(witness)
        return tuple(self.target_positions[mapping[atom]] for atom in self.source_atoms)

    def aligned_target(self, images):
        aligned = (self.target if tuple(images) == tuple(range(self.size))
                   else self.target[np.ix_(images, images)])
        return self.source if np.array_equal(aligned, self.source) else aligned

    def count(self, images):
        aligned = self.aligned_target(images)
        return int(np.count_nonzero(np.triu(np.abs(self.source - aligned) > self.thresholds, 1)))

    def incidence_lower_bound(self, seconds):
        """Relax neighbor assignments, retaining injective atom occupancy.

        Each complete mapping pays every bond event at both endpoints. The
        minimum relaxed sum of incident costs, divided by two and rounded up,
        therefore cannot exceed the actual event count. Elements and explicit
        anchors constrain the atom assignment; zero-weight contacts are included.
        """
        from mappa._group_ops import incident_event_bound
        colors = {element: i for i, element in enumerate(dict.fromkeys(self.source_elements))}
        return incident_event_bound(self.source.tolist(), self.target.tolist(), self.thresholds.tolist(),
            [colors[e] for e in self.source_elements], [colors[e] for e in self.target_elements],
            self.anchors, seconds)

    def lower_bound(self):
        """Relax atom consistency, matching pair weights within element-pair types.

        Include nonedges and subthreshold contacts. Every permitted atom mapping
        induces such a pair matching, so this relaxation cannot exceed its cost.
        Sorted interval matching maximizes the number of pairs without an event.
        """
        source, target, limits = defaultdict(list), defaultdict(list), {}
        for i, j in combinations(range(self.size), 2):
            left = tuple(sorted((self.source_elements[i], self.source_elements[j])))
            right = tuple(sorted((self.target_elements[i], self.target_elements[j])))
            source[left].append(float(self.source[i, j]))
            target[right].append(float(self.target[i, j]))
            limits[left] = float(self.thresholds[i, j])
        if {k: len(v) for k, v in source.items()} != {k: len(v) for k, v in target.items()}:
            raise ValueError('event endpoints must preserve fragment element multiplicities')
        matched = 0
        for kind, values in source.items():
            left, right = sorted(values), sorted(target[kind])
            i = j = 0
            while i < len(left) and j < len(right):
                if abs(left[i] - right[j]) <= limits[kind]:
                    matched += 1
                    i += 1
                    j += 1
                elif left[i] < right[j]:
                    i += 1
                else:
                    j += 1
        return self.size * (self.size - 1) // 2 - matched
