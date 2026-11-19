"""Compile saved fragment relationships for membership and assignment queries."""
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
    physical mapping. Partial matching and external analytics share this compiler;
    ordinary event decoding uses its prepared local responses instead.
    """
    from mappa.mappa_core.core.symmetry.constraints import FiniteConstraints
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
    program = []
    for edge_id in reversed(path.transitions):
        edge = path.graph.transitions[edge_id]
        if edge.match is not None:
            program.append(('relation', edge_id, edge.match))
        if edge.placement_generators:
            program.append(('group', edge_id, edge.placement_generators))
    values, actions = compile_actions(encoder, values, program)
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


def compile_actions(encoder, values, program):
    """Bind requested images to complete, ordered fragment and placement choices."""
    from mappa.mappa_core.core.symmetry.constraints import compile_relaxed_relation
    actions = []
    for kind, label, data in program:
        if kind == 'relation':
            relation = compile_relaxed_relation(encoder, data.relation, label, images_only=True)
            transport = {p: relation.images[r] for r, p in data.witness}
            values = [encoder.lookup(value, transport) for value in values]
            actions.append((kind, label, relation))
        else:
            start = len(encoder.program)
            values = encoder.act(values, data, ('placement', label))
            actions.append((kind, label, tuple(encoder.program[start:])))
    return values, actions


def realize_actions(encoder, actions, model, degree):
    """Return checked full permutations, retaining each local discovery certificate."""
    def number(value):
        expression = encoder.expression(value)
        return expression if isinstance(expression, int) else model.eval(
            expression, model_completion=True).as_long()

    result = []
    for kind, label, action in actions:
        if kind == 'relation':
            selected = action.realize(model)
            transport = {p: selected[r] for r, p in action.relation.witness}
            result.append(dict(kind=kind, transition=label,
                permutation=sorted(transport.items()), **action.choices(model)))
        else:
            permutation = tuple(range(degree))
            for _, _, choice, factor in action:
                selected = factor[number(choice)]
                permutation = tuple(selected[p] if p < len(selected) else p for p in permutation)
            result.append(dict(kind=kind, transition=label,
                               permutation=list(enumerate(permutation))))
    return result


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


def _tuple_orbit(seeds, moves, limit):
    """Close atom tuples under generator moves; unlisted atoms stay fixed."""
    from mappa.mappa_core.core.search.branching.construction import BranchLimitExceeded
    seen = set(seeds)
    stack = list(seen)
    while stack:
        current = stack.pop()
        for move in moves:
            image = tuple(move.get(atom, atom) for atom in current)
            if image not in seen:
                seen.add(image)
                stack.append(image)
                if len(seen) > int(limit):
                    raise BranchLimitExceeded(limit, branch_count=len(seen),
                                              stage='core_symmetry_expansion')
    return seen


def _self_moves(symmetry):
    """A fragment's own shuffles as maps between original endpoint indices."""
    atoms = symmetry.atoms
    return tuple({atom: atoms[image] for atom, image in zip(atoms, generator)
                  if atoms[image] != atom} for generator in symmetry.generators)


def _relation_images(relation, atoms, limit):
    """Images of a source tuple under the fragment's own shuffles.

    The saved witness is carried through the source fragment's shuffles and
    then the target fragment's: target shuffle o witness o source shuffle.
    Atoms move together as one tuple, so their correlations are kept. Matches
    that no shuffle of the saved witness reaches are not searched for.
    """
    witness = dict(relation.witness)
    sources = _tuple_orbit({tuple(atoms)}, _self_moves(relation.source), limit)
    placed = {tuple(witness[atom] for atom in source) for source in sources}
    return _tuple_orbit(placed, _self_moves(relation.target), limit)


def project_path_assignments(path, problem, source_atoms, *, assignment_limit):
    """Enumerate only distinct requested tuples from the exact path relation.

    The path already stores each fragment's shuffles and placements, so the
    requested atoms are carried through them directly as tuples. Choices within
    one fragment stay joint; no full mapping is expanded and no solver is used.
    """
    from mappa.mappa_core.core.search.branching.construction import BranchLimitExceeded

    source_atoms = tuple(sorted(set(map(int, source_atoms))))
    representative = path.mapping
    if any(atom not in representative for atom in source_atoms):
        return []
    if len(set(representative.values())) != len(representative):
        return []
    target_count = problem.target_atom_count
    tuples = {tuple(representative[atom] for atom in source_atoms)}
    for edge_id in reversed(path.transitions):
        edge = path.graph.transitions[edge_id]
        if edge.match is not None:
            relation = edge.match.relation
            back = {p: r for r, p in edge.match.witness}
            cache, moved = {}, set()
            for current in tuples:
                needed = tuple(dict.fromkeys(back[p] for p in current if p in back))
                if not needed:
                    moved.add(current)
                    continue
                if needed not in cache:
                    cache[needed] = _relation_images(relation, needed, assignment_limit)
                for images in cache[needed]:
                    image = dict(zip(needed, images))
                    moved.add(tuple(image[back[p]] if p in back else p for p in current))
                    if len(moved) > int(assignment_limit):
                        raise BranchLimitExceeded(assignment_limit, branch_count=len(moved),
                                                  stage='core_symmetry_expansion')
            tuples = moved
        if edge.placement_generators:
            moves = tuple({p: q for p, q in enumerate(generator) if p != q}
                          for generator in edge.placement_generators)
            tuples = _tuple_orbit(tuples, moves, assignment_limit)
    source_elements, target_elements = problem.reactant.elements, problem.product.elements
    kept = sorted(images for images in tuples if all(
        p < target_count and source_elements[r] == target_elements[p]
        for r, p in zip(source_atoms, images)))
    return [dict(zip(source_atoms, images)) for images in kept]


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
    actions = realize_actions(compiled.constraints, compiled.actions, model,
                              compiled.problem.target_atom_count)
    checked_edges = sum(len(action.relation.accepted_edges)
                        for kind, _, action in compiled.actions if kind == 'relation')
    for action in actions:
        transport = dict(action['permutation'])
        actual = {r: transport.get(p, p) for r, p in actual.items()}
    assert actual == {r: number(compiled.values[r][0]) for r in actual}
    assert len(set(actual.values())) == len(actual)
    assert all(compiled.problem.reactant.elements[r] == compiled.problem.product.elements[p]
               for r, p in actual.items())
    return dict(mapping=sorted(actual.items()), actions=actions,
                checked_fragment_edges=checked_edges, scope='full_explicit')
