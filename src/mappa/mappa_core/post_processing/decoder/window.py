"""Event windows under independent fragments and a single sweep correction."""

from dataclasses import dataclass

from mappa.mappa_core.utils.bonds import bond_event_threshold
from mappa.mappa_core.core.symmetry.event_counts import LocalEventSpectrum
from .events import FragmentEventTables


@dataclass(frozen=True)
class EventWindowSelection:
    minimum: int | None
    limit: int | None
    retained_family_ids: tuple
    retained_partitions: tuple
    lower_bounds: tuple
    branch_minima: tuple
    event_counts: tuple
    slack: int | None


@dataclass(frozen=True)
class _Composition:
    family: int
    partition: int
    responses: tuple
    boundary: int
    cut: tuple | None
    cut_owners: tuple
    local_minima: tuple
    sweep_removed: int
    sweep_width: int
    target_boundary: frozenset

    @property
    def base(self):
        return self.boundary + sum(self.local_minima)

    @property
    def lower(self):
        return self.base - self.sweep_removed


def _sums(groups, ceiling):
    """Convolve small integer supports with native big-integer bit operations."""
    values = 1
    mask = None if ceiling is None else (1 << (ceiling + 1)) - 1
    for counts in groups:
        following = 0
        for count in counts:
            following |= values << count
        values = following if mask is None else following & mask
    result = []
    while values:
        bit = values & -values
        result.append(bit.bit_length() - 1)
        values ^= bit
    return tuple(result)


def select_event_window(catalogue, index, slack, *, local):
    """Compose local counts; never optimize a complete branch mapping relation.

    The approximation counts every ordinary boundary as breaking/forming. A
    cross-fragment swept pair can remove at most two of those flat events.
    Local endpoint-conditioned spectra preserve the sweep/internal tradeoff.
    """
    compositions = []
    inventory = local.inventory
    boundary_sets = {}

    def spectrum(response, allowance, endpoint=None):
        constant = response.constant_count
        if constant is not None and endpoint is None:
            return LocalEventSpectrum(constant, (constant,), ())
        return local.counts.spectrum(
            response.mapping, response.actions,
            endpoints=() if endpoint is None else (endpoint,), allowance=allowance,
        )

    for fid, (family, partitions) in enumerate(zip(catalogue.families, inventory.uses)):
        for pid, pieces in enumerate(partitions):
            table = FragmentEventTables(local, family, pieces, inventory.cuts[fid][pid])
            responses = tuple(response for response, _ in table.fragment_responses)
            minima = tuple(
                spectrum(response, 0).minimum
                for response in responses
            )
            owners = (
                ()
                if table.sweep_cut is None
                else tuple(int(table.owners[a]) for a in table.sweep_cut)
            )
            target_boundary = table.product_boundary_set if owners else frozenset()
            target_boundary = boundary_sets.setdefault(target_boundary, target_boundary)
            compositions.append(
                _Composition(
                    fid, pid, responses, table.boundary_count, table.sweep_cut, owners, minima,
                    table.sweep_removed, table.sweep_width, target_boundary,
                )
            )
            del table
    if not compositions:
        empty = tuple(() for _ in catalogue.families)
        return EventWindowSelection(None, None, (), empty, empty, empty, empty, slack)

    paired = {}

    def sweep_counts(branch, allowance):
        # FragmentResponse contains domain dictionaries; identity is stable for
        # the owning decoder and shares equal local preparation across uses.
        key = (tuple(id(branch.responses[i]) for i in branch.cut_owners),
               branch.cut, branch.target_boundary, allowance)
        if key in paired:
            return paired[key]
        left, right = branch.cut_owners
        a, b = branch.cut
        first = spectrum(branch.responses[left], allowance, a)
        second = spectrum(branch.responses[right], allowance, b)
        tau = bond_event_threshold(
            index.elements,
            a,
            b,
            default_threshold=index.threshold,
            metal_threshold=index.metal_threshold,
        )
        source = float(index.r[a, b])
        counts = set()
        for (x,), first_counts in first.endpoint_counts:
            for (y,), second_counts in second.endpoint_counts:
                # Bins compose independent local choices. Joint occupancy is
                # validated when realizing a mapping, not inferred from a bin.
                target = float(index.p[x, y])
                counted = tuple(sorted((x, y))) in branch.target_boundary
                correction = int(source > tau) + int(counted) - int(abs(source - target) > tau)
                counts.update(i + j - correction for i in first_counts for j in second_counts)
        if not counts:
            raise ValueError("Sweep endpoint choices cannot realize the independent-fragment model")
        paired[key] = tuple(sorted(counts))
        return paired[key]

    def minimum(branch):
        if branch.cut is None:
            return branch.base
        affected = sum(branch.local_minima[i] for i in branch.cut_owners)
        return branch.base - affected + min(sweep_counts(branch, 2))

    upper = min(branch.lower + branch.sweep_width for branch in compositions)
    minima = {}
    for branch in sorted(compositions, key=lambda b: (b.lower, b.base, b.family, b.partition)):
        if branch.lower > upper:
            continue
        value = minimum(branch)
        minima[branch.family, branch.partition] = value
        upper = min(upper, value)
    limit = None if slack is None else upper + slack
    lower = [[None] * len(parts) for parts in inventory.partitions]
    minimum_rows = [[None] * len(parts) for parts in inventory.partitions]
    counts_rows = [[()] * len(parts) for parts in inventory.partitions]
    retained = [[] for _ in catalogue.families]
    for branch in compositions:
        fid, pid = branch.family, branch.partition
        lower[fid][pid] = branch.lower
        if limit is not None and branch.lower > limit:
            continue
        value = minima.get((fid, pid))
        if value is None:
            value = minimum(branch)
        minimum_rows[fid][pid] = value
        if limit is not None and value > limit:
            continue
        # Any admitted choice's excess above the independent local minima is
        # bounded by the remaining window and, at most, two cancelled events.
        allowance = None if limit is None else limit - branch.lower
        groups = [
            spectrum(response, allowance).counts
            for i, response in enumerate(branch.responses)
            if i not in branch.cut_owners
        ]
        if branch.cut is not None:
            combined = sweep_counts(branch, allowance)
            # Sweep cancellation may make its two-fragment subtotal negative.
            # Shift before bitset composition; the boundary owns those events.
            shift = min(combined)
            groups.append(tuple(count - shift for count in combined))
        else:
            shift = 0
        base = branch.boundary + shift
        counts = tuple(
            base + count for count in _sums(groups, None if limit is None else limit - base)
        )
        if counts:
            retained[fid].append(pid)
            counts_rows[fid][pid] = counts
    return EventWindowSelection(
        upper,
        limit,
        tuple(fid for fid, parts in enumerate(retained) if parts),
        tuple(map(tuple, retained)),
        tuple(map(tuple, lower)),
        tuple(map(tuple, minimum_rows)),
        tuple(map(tuple, counts_rows)),
        slack,
    )
