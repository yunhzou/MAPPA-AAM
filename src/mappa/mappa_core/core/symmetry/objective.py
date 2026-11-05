"""Original-weight event counts and a connectivity-free lower bound."""
from collections import defaultdict

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
            self.source = left.wbo[np.ix_(source.atoms, source.atoms)]
            self.target = right.wbo[np.ix_(target.atoms, target.atoms)]
            self.source_elements = tuple(left.elements[a] for a in source.atoms)
            self.target_elements = tuple(right.elements[a] for a in target.atoms)
        self.size = len(source.atoms)
        self.pairs = tuple(zip(*np.triu_indices(self.size, 1)))
        self.thresholds = np.full((self.size, self.size), float(threshold))
        if metal_threshold is not None:
            metal = np.asarray([is_metal_element(element) for element in self.source_elements], dtype=bool)
            self.thresholds[metal[:, None] | metal[None, :]] = float(metal_threshold)
        self.target_positions = {atom: i for i, atom in enumerate(target.atoms)}
        self.source_atoms = source.atoms

    def images(self, witness):
        mapping = dict(witness)
        return tuple(self.target_positions[mapping[atom]] for atom in self.source_atoms)

    def count(self, images):
        aligned = self.target[np.ix_(images, images)]
        return int(np.count_nonzero(np.triu(np.abs(self.source - aligned) > self.thresholds, 1)))

    def lower_bound(self):
        """Relax atom consistency, matching pair weights within element-pair types.

        Include nonedges and subthreshold contacts. Every permitted atom mapping
        induces such a pair matching, so this relaxation cannot exceed its cost.
        Sorted interval matching maximizes the number of pairs without an event.
        """
        source, target, limits = defaultdict(list), defaultdict(list), {}
        for i, j in self.pairs:
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
        return len(self.pairs) - matched
