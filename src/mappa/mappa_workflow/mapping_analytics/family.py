"""Reference membership in the complete saved fragment relationships."""
from functools import lru_cache, cached_property
from dataclasses import dataclass
import time

from mappa.mappa_core.core.symmetry.actions import ActionIndex


@lru_cache(maxsize=1)
def _reference_actions(degree):
    """Reuse the core's bounded immutable action domains in one endpoint frame."""
    return ActionIndex(degree)


@lru_cache(maxsize=256)
def atom_sets(degree, generators):
    """Possible atom sets only; this never asserts joint assignment compatibility."""
    parent=list(range(degree))
    def root(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]];a=parent[a]
        return a
    for generator in generators:
        for a,b in enumerate(generator):parent[root(a)]=root(b)
    groups={}
    for a in range(degree):groups.setdefault(root(a),set()).add(a)
    return tuple(frozenset(groups[root(a)]) for a in range(degree))


def _compatible_reference(path, representative, reference, source_atoms,
                          source_generators, target_generators, complete_reference,
                          source_count, target_count):
    """Reject only references outside the full ordered program's atom bounds.

    The core bounds include source self, discovery and target self choices.
    They deliberately ignore joint restrictions: surviving references still
    require the exact compiler and solver below.
    """
    requested = tuple(atom for atom in source_atoms if complete_reference or atom in reference)
    if not requested:
        return True
    index = _reference_actions(target_count)
    program = []
    for edge_id in reversed(path.transitions):
        edge = path.graph.transitions[edge_id]
        if edge.match is not None:
            program.append(index.action_domains(('relation', edge.match)))
        if edge.placement_generators:
            program.append(index.action_domains(('group', edge.placement_generators)))
    source = atom_sets(source_count, tuple(source_generators))
    target = atom_sets(target_count, tuple(target_generators))
    for atom in requested:
        wanted = reference.get(atom, target_count)
        allowed = target[wanted] if 0 <= wanted < target_count else frozenset((wanted,))
        domain = frozenset(representative.get(r, target_count) for r in source[atom])
        for action in program:
            domain = frozenset(q for p in domain
                for q in (action[p] if p < target_count else (p,)))
        if domain.isdisjoint(allowed):
            return False
    return True


def compile_path(path, problem, reference, *, source_atoms, source_generators=(),
                 target_generators=(), complete_reference=True):
    """Compile saved local relations, placement actions and the reference.

    Every local discovery and both self rearrangements remain joint choices.
    Endpoint equivalence normalizes reference identity only, never the returned
    physical mapping. Decoding does not call this reference query.
    """
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints, compile_relaxed_relation
    started = time.perf_counter()
    encoder = FiniteConstraints()
    source_atoms = tuple(source_atoms)
    nr, np_ = problem.source_atom_count, problem.target_atom_count
    representative = path.mapping
    values = [encoder.constant(representative.get(r, np_)) for r in range(nr)]
    actions = []
    if not _compatible_reference(path, representative, reference, source_atoms,
                                 source_generators, target_generators,
                                 complete_reference, nr, np_):
        encoder.failed = True
        return CompiledFamily(encoder, values, representative, actions, path, problem,
                              time.perf_counter() - started)
    for edge_id in reversed(path.transitions):
        edge = path.graph.transitions[edge_id]
        fragment = edge.match
        if fragment is not None:
            relation = compile_relaxed_relation(encoder, fragment.relation, edge_id)
            transport = {p: relation.images[r] for r, p in fragment.witness}
            values = [encoder.lookup(value, transport) for value in values]
            actions.append(('relation', edge_id, relation))
        if edge.placement_generators:
            start = len(encoder.program)
            values = encoder.act(values, edge.placement_generators, ('placement', edge_id))
            actions.append(('group', edge_id, tuple(encoder.program[start:])))
    if len(set(representative.values())) != len(representative):
        encoder.failed = True
    for r in representative:
        encoder.restrict(values[r][0], {p for p in values[r][1]
            if p < np_ and problem.reactant.elements[r] == problem.product.elements[p]})
    physical_values = values
    source_choices = encoder.act([encoder.constant(r) for r in source_atoms],
                                 source_generators, 'source_equivalence')
    normalized = [encoder.lookup(value, dict(enumerate(physical_values))) for value in source_choices]
    normalized = encoder.act(normalized, target_generators, 'target_equivalence')
    for atom, value in zip(source_atoms, normalized):
        if complete_reference or atom in reference:
            encoder.restrict(value[0], {reference.get(atom, np_)})
    encoder.propagate()
    return CompiledFamily(encoder, physical_values, representative, actions, path, problem,
                          time.perf_counter() - started)


@dataclass
class CompiledFamily:
    constraints: object
    physical_values: list
    representative: dict
    actions: list
    path: object
    problem: object
    encoding_seconds: float

    @cached_property
    def values(self):
        """Build output expressions only when a caller needs physical images."""
        return [(self.constraints.expression(value), frozenset(self.constraints.domain(value)))
                for value, _ in self.physical_values]

    @cached_property
    def solver(self):
        """Materialize a joint query only for callers adding further constraints."""
        import z3
        solver=z3.Solver()
        self.constraints.emit(solver)
        return solver

    def realize(self, model):
        return _realize(self,model)


def project_path_assignments(path, problem, source_atoms, *, assignment_limit):
    """Enumerate only distinct requested tuples from the exact path relation."""
    import z3
    from mappa.mappa_core.core.search.branching.construction import BranchLimitExceeded

    source_atoms = tuple(sorted(set(map(int, source_atoms))))
    if any(atom not in path.mapping for atom in source_atoms):
        return []
    compiled = compile_path(path, problem, {}, source_atoms=source_atoms,
                            complete_reference=False)
    values = tuple(compiled.values[atom][0] for atom in source_atoms)
    assignments = []
    while True:
        status = compiled.solver.check()
        if status == z3.unsat:
            break
        if status == z3.unknown:
            raise RuntimeError(
                f"Core assignment query was inconclusive: {compiled.solver.reason_unknown()}")
        model = compiled.solver.model()
        images = tuple(value if isinstance(value, int)
                       else model.eval(value, model_completion=True).as_long()
                       for value in values)
        assignments.append(dict(zip(source_atoms, images)))
        if len(assignments) > int(assignment_limit):
            raise BranchLimitExceeded(assignment_limit, branch_count=len(assignments),
                                      stage='core_symmetry_expansion')
        compiled.solver.add(z3.Or(*(value != image
            for value, image in zip(values, images))))
    return sorted(assignments, key=lambda mapping: tuple(mapping.values()))


def query_path(path, problem, reference, *, source_atoms, source_generators=(),
               target_generators=(), complete_reference=True, timeout_ms=10000):
    """Return a certified full explicit-atom realization, UNSAT, or unknown."""
    import z3
    start, cpu = time.perf_counter(), time.process_time()
    compiled=compile_path(path,problem,reference,source_atoms=source_atoms,
        source_generators=source_generators,target_generators=target_generators,
        complete_reference=complete_reference)
    encoding = dict(wall_seconds=time.perf_counter()-start, cpu_seconds=time.process_time()-cpu)
    left=timeout_ms-1000*encoding['wall_seconds']
    if left<=0:
        return 'unknown',dict(reason='encoding budget',seconds=encoding['wall_seconds'],
                              timing=dict(encoding=encoding))
    solve, solve_cpu = time.perf_counter(), time.process_time()
    status, model, checks=compiled.constraints.solve(max(1,int(left)))
    timing = dict(encoding=encoding, solving=dict(wall_seconds=time.perf_counter()-solve,
                                                 cpu_seconds=time.process_time()-solve_cpu))
    if status!=z3.sat:
        return ('not_recovered' if status==z3.unsat else 'unknown'),dict(
            reason='unsat' if status==z3.unsat else 'reference component budget',
            seconds=time.perf_counter()-start,timing=timing,checks=checks)
    result=compiled.realize(model)
    result['seconds']=time.perf_counter()-start
    result['timing']=timing
    result['checks']=checks
    return 'recovered',result


def _realize(compiled, model):
    def number(value):
        return value if isinstance(value, int) else model.eval(value, model_completion=True).as_long()

    actual = dict(compiled.representative)
    actions = []
    checked_edges = 0
    for kind, edge_id, action in compiled.actions:
        if kind == 'relation':
            selected = action.realize(model)
            transport = {p: selected[r] for r, p in action.relation.witness}
            actions.append(dict(kind=kind, transition=edge_id,
                permutation=sorted(transport.items()), **action.choices(model)))
            checked_edges += len(action.relation.accepted_edges)
        else:
            permutation = tuple(range(compiled.problem.target_atom_count))
            for _, _, choice, factor in action:
                selected = factor[number(compiled.constraints.expression(choice))]
                permutation = tuple(selected[p] if p < len(selected) else p for p in permutation)
            transport = dict(enumerate(permutation))
            actions.append(dict(kind=kind, transition=edge_id,
                                permutation=sorted(transport.items())))
        actual = {r: transport.get(p, p) for r, p in actual.items()}
    assert actual == {r: number(compiled.values[r][0]) for r in actual}
    assert len(set(actual.values())) == len(actual)
    assert all(compiled.problem.reactant.elements[r] == compiled.problem.product.elements[p]
               for r, p in actual.items())
    return dict(mapping=sorted(actual.items()), actions=actions,
                checked_fragment_edges=checked_edges, scope='full_explicit')
