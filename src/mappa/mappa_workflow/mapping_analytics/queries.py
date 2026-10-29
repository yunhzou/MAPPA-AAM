"""Explicit assignment queries over already decoded mapping families."""
from dataclasses import dataclass
import math
import time

@dataclass(frozen=True)
class ShuffleQuery:
    status: str  # allowed / forbidden / unknown
    mapping: dict | None = None
    actions: tuple = ()
    family_id: int | None = None
    reason: str = ""


def query_assignments(decoded, candidate, assignments, *, same_events=True, seconds=10.0):
    """Can these R->P assignments coexist in a retained mapping family?

    Unspecified atoms may move jointly. To test exactly one swap, supply
    the entire modified witness. By default retain the candidate's event
    class. Search the converted catalogue, never the old free-choice program. Negative answers concern this saved catalogue only. Budgets
    are soft (encoding and native calls are not preemptible).
    """
    import z3
    from mappa.mappa_workflow.mapping_analytics.family import compile_path
    from mappa.mappa_workflow.mapping_analytics.event_query import _event_model, _same_event_orbit

    decoded._check_candidate(candidate)
    if not math.isfinite(seconds) or seconds < 0:
        raise ValueError("seconds must be finite and nonnegative")
    fixed = dict(assignments)
    problem = decoded.catalogue.problem
    if any(
        not isinstance(a, int)
        or not isinstance(b, int)
        or a not in range(problem.source_atom_count)
        or b not in range(problem.target_atom_count)
        for a, b in fixed.items()
    ):
        raise ValueError("assignments require valid integer endpoint indices")
    if len(set(fixed.values())) != len(fixed) or any(
        problem.reactant.elements[a] != problem.product.elements[b]
        for a, b in fixed.items()
    ):
        return ShuffleQuery(
            "forbidden", reason="noninjective or element-incompatible assignments"
        )
    deadline = time.perf_counter() + seconds
    pattern = decoded.index.describe(
        decoded.index.vector(candidate.mapping)
    )
    unknown = False
    for i, family in enumerate(decoded.catalogue.families):
        if time.perf_counter() >= deadline:
            return ShuffleQuery("unknown", reason="time budget")
        compiled = compile_path(
            family.as_path(problem),
            problem,
            fixed,
            source_atoms=tuple(fixed),
            complete_reference=False,
        )
        if same_events:
            try:
                terms, total, _ = _event_model(compiled, decoded.index, deadline)
            except TimeoutError:
                return ShuffleQuery("unknown", reason="event encoding budget")
            compiled.solver.add(
                _same_event_orbit(terms, total, pattern, decoded.index)
            )
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            return ShuffleQuery("unknown", reason="encoding budget")
        compiled.solver.set(timeout=max(1, int(remaining * 1000)))
        status = compiled.solver.check()
        if status == z3.sat:
            witness = compiled.realize(compiled.solver.model())
            mapping = dict(witness["mapping"])
            assert all(mapping[a] == b for a, b in fixed.items())
            if same_events:
                assert (
                    decoded.index.describe(decoded.index.vector(mapping))[
                        "id"
                    ]
                    == candidate.id
                )
            return ShuffleQuery("allowed", mapping, tuple(witness["actions"]), i)
        unknown |= status == z3.unknown
    return ShuffleQuery(
        "unknown" if unknown else "forbidden", reason="saved-family query"
    )
