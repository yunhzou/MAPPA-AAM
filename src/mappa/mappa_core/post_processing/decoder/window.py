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
        best = None
        for lower, fid, partition, cut in sorted(candidates):
            if best is not None and lower >= best:
                break
            table = self._table(decoded, fid, partition, cut)
            best = table.count_bounds[1] if best is None else min(best, table.count_bounds[1])
            if table.count_bounds[0] < best:
                best = min(best, self._minimum(decoded, table, best))
            del table
        self.minimum = best
        self.maximum = None if best is None else best + decoded.event_slack
        self.contexts = []
        self.bounds = {}
        for lower, fid, partition, cut in candidates:
            if lower > self.maximum:
                continue
            table = self._table(decoded, fid, partition, cut)
            # Keep conditional local responses inside the window. A bound need
            # only reject impossible compositions; it must not turn export into
            # an eager minimum/support query for every retained fragment.
            lower = table.count_bounds[0]
            if lower <= self.maximum:
                key = fid, partition, cut
                self.contexts.append(key)
                self.bounds[key] = lower
            del table
        self.contexts = tuple(self.contexts)

    @staticmethod
    def _table(decoded, fid, partition, cut):
        key = decoded._table_key(fid, partition, cut)
        return decoded._prepare_table(fid, partition, key)

    @staticmethod
    def _minimum(decoded, table, upper):
        boundary = table.boundary_events
        responses = tuple(response for response, _ in table.fragment_responses)
        minima = tuple(decoded.local.minimum(response) for response in responses)
        lower = boundary.lower + sum(minima)
        if lower >= upper or table.sweep_cut is None:
            return lower
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
        # Most swept pairs have a fixed response even when their endpoint
        # fragments contain many choices. Certify that before a joint query.
        lower += reserve + min(cost for (p, q), cost in costs.items() if p != q)
        if lower >= upper:
            return lower
        from .local_counts import compile_local
        query = compile_local(decoded.index, selected)
        x, y = query.images[0][a], query.images[1][b]
        query.solver.add(x != y)
        correction = query.pair_cost(x, y, costs)
        base = boundary.lower + sum(value for i, value in enumerate(minima) if i not in (left, right))
        expression = base + reserve + sum(query.counts) + correction
        return query.minimum(expression, lower, table.count_bounds[1])
