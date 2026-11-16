"""Bound branch preparation by attained witnesses; retain incomplete minima."""
from dataclasses import replace

from .events import prepare_boundary_events


def _choice_key(value):
    """Intern equal choices while retaining their original action-local meaning."""
    if isinstance(value, dict):
        return 'mapping', tuple(sorted((key, _choice_key(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return 'sequence', tuple(_choice_key(item) for item in value)
    return value


def _restore_choice(choice):
    if isinstance(choice, dict):
        return {key: ({int(a): int(b) for a, b in value.items()} if key == 'discovery'
                      else tuple(value)) for key, value in choice.items()}
    return None if choice is None else tuple(choice)


class EventWindow:
    """Window of the retained binary-symmetry, independent-boundary model.

    Local witnesses establish an attained upper bound. When their minima compose,
    it lies within the smallest pre-sweep count plus two. An unfinished local
    search instead keeps its honest lower bound and the attained upper bound.
    """
    def __init__(self, decoded):
        candidates = []
        for fid, partitions in enumerate(decoded.cut_contexts):
            for partition, cuts in enumerate(partitions):
                for cut, bonds in enumerate(cuts):
                    boundary = prepare_boundary_events(decoded.index, decoded.catalogue.families[fid].mapping,
                        decoded.inventory.uses[fid][partition], bonds)
                    candidates.append((boundary.lower, (fid, partition, cut)))
        bounds, witnesses, responses = {}, {}, {}
        upper = None
        for cheap, key in sorted(candidates):
            if upper is not None and cheap > upper + decoded.event_slack:
                bounds[key] = cheap
                continue
            table = self._table(decoded, *key)
            responses[key] = tuple(response for response, _ in table.fragment_responses)
            bound = table.boundary_events.lower
            for response, _ in table.fragment_responses:
                counts = response.event_counts
                bound += counts.minimum if counts.minimum_proven else counts.lower
            bounds[key] = bound
            witness = table.composed_minimum()
            if witness is None:
                import numpy as np
                program = table.program
                prefix = np.arange(program.size)
                count = program.fixed_events + int(program.score(prefix, program.pairs))
                mapping, choices = program.query_result(prefix, count, 0., None)
                witness = count, mapping, choices
            assert bound <= witness[0], 'a local lower bound exceeded an attained branch count'
            witnesses[key] = witness
            upper = witness[0] if upper is None else min(upper, witness[0])
            del table
        self.initial_lower = min(bounds.values(), default=None)
        self.initial_maximum = upper
        self.results = {}
        best = upper
        for key in sorted(witnesses, key=lambda key: (bounds[key], key)):
            if bounds[key] >= best:
                continue
            table = self._table(decoded, *key)
            answer = table.minimum_result
            known, mapping, choices = witnesses[key]
            if known < answer.event_count:
                answer = replace(answer, event_count=known, witness=mapping, choices=choices,
                                 complete=known == bounds[key])
            self.results[key] = answer
            bounds[key] = answer.event_count if answer.complete else max(bounds[key], answer.lower_bound)
            best = min(best, answer.event_count)
            del table
        self.minimum = best
        self.lower_bound = min(bounds.values(), default=None)
        self.minimum_proven = best is not None and self.lower_bound == best
        self.maximum = None if best is None else best + decoded.event_slack
        self.contexts = tuple(key for _, key in candidates if bounds[key] <= self.maximum)
        self.bounds = {key: bounds[key] for key in self.contexts}
        # The witnesses already exist. Persist shared action choices, not a full
        # atom mapping for every branch, and do not launch additional searches.
        self.branch_minima, self.minimum_choices = {}, []
        identities, values = {}, {}
        for key in self.contexts:
            answer = self.results.get(key)
            count, _, choices = witnesses[key]
            if answer is not None:
                count, choices = answer.event_count, answer.choices
            references = []
            for choice in choices:
                identity = id(choice)
                if identity not in identities:
                    frozen = _choice_key(choice)
                    if frozen not in values:
                        values[frozen] = len(self.minimum_choices)
                        self.minimum_choices.append(choice)
                    identities[identity] = values[frozen]
                references.append(identities[identity])
            self.branch_minima[key] = int(count), int(bounds[key]), tuple(references)
        self.minimum_choices = tuple(self.minimum_choices)
        retained = {id(response) for key in self.contexts for response in responses[key]}
        decoded.local.retain_responses(retained)

    @staticmethod
    def from_record(record):
        window = EventWindow.__new__(EventWindow)
        for name in ('minimum', 'maximum', 'lower_bound', 'minimum_proven',
                     'initial_lower', 'initial_maximum'):
            setattr(window, name, record[name])
        window.bounds = {tuple(key): value for key, value in record['bounds']}
        window.contexts = tuple(tuple(key) for key in record['contexts'])
        window.results = {}
        # These additive fields are absent from previously published v14 files.
        # Such files retain their existing query behavior; no minima are invented.
        window.branch_minima = {tuple(key): (int(count), int(lower), tuple(choices))
                               for key, (count, lower, choices) in record.get('branch_minima', ())}
        window.minimum_choices = tuple(_restore_choice(choice)
                                       for choice in record.get('minimum_choices', ()))
        return window

    def to_record(self):
        return {**{name: getattr(self, name) for name in (
            'minimum', 'maximum', 'lower_bound', 'minimum_proven', 'initial_lower', 'initial_maximum')},
            'bounds': tuple(self.bounds.items()), 'contexts': self.contexts,
            'branch_minima': tuple(self.branch_minima.items()),
            'minimum_choices': self.minimum_choices}

    @staticmethod
    def _table(decoded, fid, partition, cut):
        key = decoded._table_key(fid, partition, cut)
        return decoded._prepare_table(fid, partition, key)
