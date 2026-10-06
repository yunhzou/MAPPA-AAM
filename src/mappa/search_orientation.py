"""Direction selection outside the directed, symmetry-compressed AAM engine."""
from dataclasses import dataclass, replace

from .domain import AAMProblem, AAMSearchConfig, AAMResult


@dataclass(frozen=True)
class AAMSearchPlan:
    """Keep input R/P identities separate from search source/target identities.

    The raw AAM result and its symmetry actions retain search orientation.
    Convert only a concrete mapping after selecting/realizing its search path.
    Do not invert generators or reinterpret a compressed graph in place.
    """
    input_problem: AAMProblem
    problem: AAMProblem
    config: AAMSearchConfig
    reversed: bool

    def to_input_mapping(self, mapping):
        pairs = dict(mapping)
        return {r:p for p,r in pairs.items()} if self.reversed else pairs

    def to_search_mapping(self, mapping):
        return self.to_input_mapping(mapping)

    @property
    def direction(self):
        return 'P_to_R' if self.reversed else 'R_to_P'


def plan_aam_search(problem, config=None):
    """Search smaller explicit-atom endpoint first; retain R→P on ties."""
    config = config or AAMSearchConfig()
    reverse = problem.source_atom_count > problem.target_atom_count
    search_problem = AAMProblem(problem.product, problem.reactant, problem.name) if reverse else problem
    search_config = replace(config, anchors=tuple((p,r) for r,p in config.anchors)) if reverse else config
    return AAMSearchPlan(problem, search_problem, search_config, reverse)


@dataclass(frozen=True)
class DirectedAAMResult:
    """A directed raw search plus conversion back to the user's atom indices."""
    plan: AAMSearchPlan
    aam: AAMResult

    def to_input_mapping(self, mapping):
        return self.plan.to_input_mapping(mapping)


def search_aam_directions(problem, config=None, *, direction='both', workers=1,
                          intermediate_dir=None, **search_options):
    """Run forward, reverse, smaller-first, larger-first, or both searches.

    Results retain native directed families; only concrete witnesses are
    inverted. Both directions run sequentially and share the worker ceiling.
    Search anchors are supplied in input R->P indices and inverted as needed.
    Post-processing and cross-direction output identity remain separate.
    """
    from pathlib import Path
    from .aam import search_aam
    config = config or AAMSearchConfig()
    choices = {'forward': (False,), 'reverse': (True,), 'both': (False, True),
               'smaller_first': (problem.source_atom_count > problem.target_atom_count,),
               'larger_first': (problem.source_atom_count <= problem.target_atom_count,)}
    if direction not in choices:
        raise ValueError(f'unknown direction: {direction}')
    results = []
    for reverse in choices[direction]:
        directed = AAMProblem(problem.product, problem.reactant, problem.name) if reverse else problem
        cfg = replace(config, anchors=tuple((p,r) for r,p in config.anchors)) if reverse else config
        plan = AAMSearchPlan(problem, directed, cfg, reverse)
        folder = None if intermediate_dir is None else Path(intermediate_dir) / plan.direction
        result = search_aam(directed, cfg, workers=workers, intermediate_dir=folder, **search_options)
        results.append(DirectedAAMResult(plan, result))
    return tuple(results)
