"""Select additive fragment compositions before resolving their sweep bond."""
from .events import prepare_boundary_events


class EventWindow:
    """Minimum and inclusive window of the independent-fragment event model."""
    def __init__(self, decoded):
        candidates = []
        for fid, partitions in enumerate(decoded.cut_contexts):
            for partition, cuts in enumerate(partitions):
                for cut, bonds in enumerate(cuts):
                    boundary = prepare_boundary_events(decoded.index, decoded.catalogue.families[fid].mapping,
                        decoded.inventory.uses[fid][partition], bonds)
                    candidates.append((boundary.lower, fid, partition, cut))
        ordered = sorted(candidates)
        bounds = {}
        lower_minimum = None
        unresolved = []
        witness_minimum = None
        for lower, fid, partition, cut in ordered:
            if lower_minimum is not None and lower >= lower_minimum:
                break
            table = self._table(decoded, fid, partition, cut)
            if lower_minimum is not None and table.count_bounds[0] >= lower_minimum:
                del table
                continue
            lower = table.count_bounds[0]
            bounds[fid, partition, cut] = lower
            witness_minimum = (table.count_bounds[1] if witness_minimum is None else
                               min(witness_minimum, table.count_bounds[1]))
            if all(response.event_counts.lower == response.event_counts.bounds[1]
                   for response, _ in table.fragment_responses):
                lower_minimum = lower if lower_minimum is None else min(lower_minimum, lower)
            else:
                unresolved.append((lower, fid, partition, cut))
            del table
        for lower, fid, partition, cut in sorted(unresolved):
            if lower_minimum is not None and lower >= lower_minimum:
                break
            table = self._table(decoded, fid, partition, cut)
            lower = table.boundary_events.lower + sum(self._local_minima(decoded, table))
            bounds[fid, partition, cut] = lower
            lower_minimum = lower if lower_minimum is None else min(lower_minimum, lower)
            del table
        initial_maximum = None if lower_minimum is None else lower_minimum + 2
        best = (witness_minimum if witness_minimum is not None and
                witness_minimum <= initial_maximum else None)
        for lower, fid, partition, cut in ordered:
            if lower > initial_maximum or (best is not None and lower >= best):
                break
            table = self._table(decoded, fid, partition, cut)
            upper = initial_maximum if best is None else min(initial_maximum, best - 1)
            if bounds.get((fid, partition, cut), table.count_bounds[0]) > upper:
                del table
                continue
            minima = self._local_minima(decoded, table)
            lower = table.boundary_events.lower + sum(minima)
            bounds[fid, partition, cut] = lower
            if lower <= upper:
                value = (lower if table.sweep_cut is None else
                         self._sweep_minimum(decoded, table, minima, upper))
                if value is not None:
                    best = value if best is None else min(best, value)
            del table
            if best == lower_minimum:
                break
        if candidates and best is None:
            raise ValueError('No attainable event count in the initial lower-bound + 2 window')
        self.minimum = best
        self.maximum = None if best is None else best + decoded.event_slack
        self.contexts = []
        self.bounds = {}
        for lower, fid, partition, cut in candidates:
            if lower > self.maximum:
                continue
            key = fid, partition, cut
            if key not in bounds:
                table = self._table(decoded, fid, partition, cut)
                # Reporting keeps conditional responses, so it needs a safe
                # exclusion bound rather than another minimum query.
                bounds[key] = table.count_bounds[0]
                del table
            lower = bounds[key]
            if lower <= self.maximum:
                self.contexts.append(key)
                self.bounds[key] = lower
        self.contexts = tuple(self.contexts)

    @staticmethod
    def _table(decoded, fid, partition, cut):
        key = decoded._table_key(fid, partition, cut)
        return decoded._prepare_table(fid, partition, key)

    @staticmethod
    def _local_minima(decoded, table):
        return tuple(response.event_counts.minimum for response, _ in table.fragment_responses)

    @staticmethod
    def _sweep_minimum(decoded, table, minima, upper):
        """Resolve only the sweep endpoints inside the initial count window."""
        boundary = table.boundary_events
        responses = tuple(response for response, _ in table.fragment_responses)
        lower = boundary.lower + sum(minima)
        a, b = table.sweep_cut
        left, right = int(table.owners[a]), int(table.owners[b])
        selected = responses[left], responses[right]
        term = next(iter(table.boundary.values()))[0]
        costs = {(p, q): (int(abs(float(decoded.index.p[p, q]) - term.response.source_weight)
                                 > term.response.threshold)
                 - int(tuple(sorted((p, q))) in table.product_boundary_set)) if p != q else 0
                 for p in selected[0].atom_domains[a] for q in selected[1].atom_domains[b]
                 }
        reserve = boundary.sweep_removed - int(boundary.source_cut_counted)
        # Fixed sweep responses need only the already prepared local counts.
        lower += reserve + min(cost for (p, q), cost in costs.items() if p != q)
        if lower > upper:
            return None
        if lower == table.count_bounds[1]:
            return lower
        base = boundary.lower + sum(value for i, value in enumerate(minima) if i not in (left, right))
        def compose(xs, ys):
            return min((base + reserve + cx + cy + costs[x, y]
                        for x, cx in xs.items() for y, cy in ys.items() if x != y), default=None)
        # A retained valid choice can certify the lower endpoint without any
        # endpoint contraction. Otherwise use the shared conditional responses.
        attained = compose(selected[0].event_counts.witness_endpoints(a),
                           selected[1].event_counts.witness_endpoints(b))
        if attained == lower:
            return attained
        lower = compose(selected[0].event_counts.endpoint_bounds(a),
                        selected[1].event_counts.endpoint_bounds(b))
        if lower is None or lower > upper:
            return None
        if attained == lower:
            return attained
        value = compose(selected[0].event_counts.endpoint(a), selected[1].event_counts.endpoint(b))
        return value if value is not None and value <= upper else None
