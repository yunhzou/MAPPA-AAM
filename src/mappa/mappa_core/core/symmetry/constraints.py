"""Finite correlated constraints, reduced before constructing the exact solver query.

Each permutation factor has one shared choice variable. Local propagation only
removes values without support; the solver retains all correlations afterward.
"""
from collections import deque, OrderedDict
from dataclasses import dataclass
from functools import lru_cache
import time


@lru_cache(maxsize=512)
def group_factors(generators):
    from .groups import permutation_group
    if not generators:return ()
    group=permutation_group(generators, len(generators[0]))
    return tuple(tuple(tuple(p) for _, p in level)
                 for level in reversed(group.transversals()) if len(level)>1)


class SymbolicActions:
    """Shared compact permutation-action encoder for feasibility and quotient queries."""
    def __init__(self,solver,*,fresh=False):
        self.solver=solver;self.fresh=fresh;self.variables=[];self.program=[]

    @staticmethod
    def constant(v):return (int(v),frozenset((int(v),)))

    def symbol(self,domain):
        import z3
        domain=tuple(sorted(set(domain)))
        if len(domain)==1:return self.constant(domain[0])
        value=z3.FreshInt('orbit') if self.fresh else z3.Int(f'q{len(self.variables)}')
        self.variables.append(value)
        self.solver.add(z3.Or(*(value==p for p in domain)))
        return value,frozenset(domain)

    def lookup(self,value,table):
        import z3
        expression,support=value
        if isinstance(expression,int):return table.get(expression,value)
        images=[(p,table.get(p,self.constant(p))) for p in sorted(support)]
        if all(isinstance(v[0],int) and v[0]==p for p,v in images):return value
        out=images[-1][1][0]
        for p,v in reversed(images[:-1]):out=z3.If(expression==p,v[0],out)
        return out,frozenset().union(*(v[1] for _,v in images))

    def act(self,values,generators,label):
        for factor in group_factors(tuple(tuple(g) for g in generators)):
            active=set().union(*(v[1] for v in values))
            moved={p for p in active if p<len(factor[0]) and any(g[p]!=p for g in factor)}
            if not moved:continue
            choice=self.symbol(range(len(factor)))
            table={p:self.lookup(choice,{i:self.constant(g[p]) for i,g in enumerate(factor)}) for p in moved}
            values=[self.lookup(v,table) for v in values]
            self.program.append(('group',label,choice[0],factor))
        return values


# Keys contain the complete reduced constraints, including the requested output.
# Only conclusive answers are reused; bounded to keep worker memory predictable.
_component_cache = OrderedDict()


@dataclass(frozen=True)
class Choice:
    index: int


class FiniteConstraints:
    def __init__(self):
        self.domains=[]
        self.constraints=[]
        self.program=[]
        self.failed=False

    @staticmethod
    def constant(value):
        return int(value), frozenset((int(value),))

    def domain(self, value):
        return self.domains[value.index] if isinstance(value,Choice) else {value}

    def symbol(self, domain):
        domain=set(domain)
        if len(domain)==1:return self.constant(next(iter(domain)))
        value=Choice(len(self.domains))
        self.domains.append(domain)
        return value,frozenset(domain)

    def restrict(self, value, allowed):
        remaining=self.domain(value).intersection(allowed)
        if not remaining:self.failed=True
        if isinstance(value,Choice) and remaining!=self.domains[value.index]:
            self.domains[value.index]=remaining
            return {value.index}
        return set()

    def relation(self, values, rows):
        self.constraints.append(('table',tuple(values),tuple(rows)))

    def distinct(self, values):
        self.constraints.append(('distinct',tuple(values),None))

    def lookup(self, value, table):
        expression,support=value
        if isinstance(expression,int):return table.get(expression,value)
        entries={p:table.get(p,self.constant(p))[0] for p in support}
        if all(v==p for p,v in entries.items()):return value
        output=self.symbol(set().union(*(self.domain(v) for v in entries.values())))
        self.constraints.append(('lookup',(expression,output[0]),entries))
        return output

    def act(self, values, generators, label):
        for factor in group_factors(tuple(tuple(g) for g in generators)):
            active=set().union(*(v[1] for v in values))
            if not any(p<len(factor[0]) and any(g[p]!=p for g in factor) for p in active):continue
            choice=self.symbol(range(len(factor)))
            values=self.apply_factor(values,factor,choice[0])
            self.program.append(('group',label,choice[0],factor))
        return values

    def apply_factor(self, values, factor, choice):
        """Apply one shared factor choice, also reused for its inverse."""
        transformed=[]
        for value in values:
            images={g[p] if p<len(g) else p for p in value[1] for g in factor}
            if all(all(p>=len(g) or g[p]==p for g in factor) for p in value[1]):
                transformed.append(value);continue
            output=self.symbol(images)
            # Restrict input, output and shared choice before constructing rows.
            self.constraints.append(('action',(value[0],choice,output[0]),factor))
            transformed.append(output)
        return transformed

    def action_rows(self, values, factor):
        source,choice,target=values
        allowed=self.domain(target)
        for p in self.domain(source):
            for c in self.domain(choice):
                g=factor[c];q=g[p] if p<len(g) else p
                if q in allowed:
                    yield p,c,q

    def propagate(self, check=None):
        neighbors=[set() for _ in self.domains]
        for n,(kind,values,data) in enumerate(self.constraints):
            for value in (*values,*(data.values() if kind=='lookup' else ())):
                if isinstance(value,Choice):neighbors[value.index].add(n)
        # Reference assignments are attached last: propagate them backwards first.
        queue=deque(reversed(range(len(self.constraints))));queued=set(queue)
        while queue and not self.failed:
            if check is not None:check()
            n=queue.popleft();queued.remove(n)
            kind,values,data=self.constraints[n];changed=set()
            if kind in ('table','action'):
                valid=[]
                for row in self.action_rows(values,data) if kind=='action' else data:
                    assignments={}
                    for value,item in zip(values,row):
                        if item not in self.domain(value) or (value in assignments and assignments[value]!=item):break
                        assignments[value]=item
                    else:valid.append(row)
                if not valid:self.failed=True;break
                for i,value in enumerate(values):changed |= self.restrict(value,{r[i] for r in valid})
                if kind=='table':self.constraints[n]=(kind,values,tuple(valid))
            elif kind=='lookup':
                selector,output=values
                eligible={p for p in self.domain(selector) if self.domain(output)&self.domain(data[p])}
                changed |= self.restrict(selector,eligible)
                changed |= self.restrict(output,set().union(*(self.domain(data[p]) for p in eligible)))
                if len(eligible)==1:
                    changed |= self.restrict(data[next(iter(eligible))],self.domain(output))
            else:
                fixed=[next(iter(self.domain(v))) for v in values if len(self.domain(v))==1]
                if len(fixed)!=len(set(fixed)):self.failed=True;break
                for value in values:
                    if len(self.domain(value))>1:changed |= self.restrict(value,self.domain(value)-set(fixed))
            for index in changed:
                for other in neighbors[index]-queued:queue.append(other);queued.add(other)

    def predicate(self, value, item):
        import z3
        domain=self.domain(value)
        if item not in domain:return False
        if len(domain)==1:return True
        return z3.Bool(f'q{value.index}_{item}')

    def expression(self, value):
        import z3
        domain=self.domain(value)
        if len(domain)==1:return next(iter(domain))
        return z3.Sum(*(z3.If(self.predicate(value,item),item,0) for item in sorted(domain)))

    def emit(self, solver):
        """Send the reduced Boolean constraints through one native parser call."""
        from itertools import product
        self.propagate()
        if self.failed:solver.add(False);return
        statements=[]

        def predicate(value, item):
            domain=self.domain(value)
            if item not in domain:return 'false'
            if len(domain)==1:return 'true'
            return f'q{value.index}_{item}'

        def negate(flag):
            return {'true':'false','false':'true'}.get(flag, f'(not {flag})')

        def clause(flags):
            if 'true' in flags:return
            flags=[flag for flag in flags if flag!='false']
            expression=(f"(or {' '.join(flags)})" if len(flags)>1
                        else flags[0] if flags else 'false')
            statements.append(f'(assert {expression})')

        for i,domain in enumerate(self.domains):
            if len(domain)>1:
                flags=[predicate(Choice(i),v) for v in sorted(domain)]
                statements.extend(f'(declare-fun {flag} () Bool)' for flag in flags)
                statements.append(f"(assert ((_ pbeq 1 {' '.join('1' for _ in flags)}) {' '.join(flags)}))")
        for kind,values,data in self.constraints:
            if kind in ('table','action'):
                # A local implication per input tuple avoids asking the solver
                # to select from a large disjunction of complete assignments.
                supported={}
                for row in self.action_rows(values,data) if kind=='action' else data:
                    supported.setdefault(row[:-1],set()).add(row[-1])
                for inputs in product(*(sorted(self.domain(v)) for v in values[:-1])):
                    premise=[negate(predicate(v,x)) for v,x in zip(values[:-1],inputs)]
                    outputs=[predicate(values[-1],x) for x in sorted(supported.get(inputs,()))]
                    clause(premise+outputs)
            elif kind=='lookup':
                selector,output=values
                for item in sorted(self.domain(selector)):
                    entry=data[item]
                    for result in sorted(self.domain(output)|self.domain(entry)):
                        left,right=predicate(output,result),predicate(entry,result)
                        if left!=right:
                            statements.append(f'(assert (=> {predicate(selector,item)} (= {left} {right})))')
            else:
                for item in sorted(set().union(*(self.domain(v) for v in values))):
                    flags=[predicate(v,item) for v in values]
                    fixed=flags.count('true')
                    choices=[flag for flag in flags if flag not in ('true','false')]
                    if fixed>1:statements.append('(assert false)')
                    elif choices:
                        statements.append(f"(assert ((_ at-most {1-fixed}) {' '.join(choices)}))")
        solver.from_string('\n'.join(statements))

    def components(self):
        """Separate only choices with no remaining shared constraint."""
        parent=list(range(len(self.domains)))
        def root(i):
            while parent[i]!=i:
                parent[i]=parent[parent[i]];i=parent[i]
            return i
        records=[]
        for kind,values,data in self.constraints:
            if kind=='lookup':
                data={p:data[p] for p in self.domain(values[0])}
            dependencies=(*values,*(data.values() if kind=='lookup' else ()))
            active=sorted({v.index for v in dependencies
                           if isinstance(v,Choice) and len(self.domain(v))>1})
            if not active:continue
            for i in active[1:]:parent[root(i)]=root(active[0])
            records.append((active,(kind,values,data)))
        groups={}
        for active,record in records:
            ids,rows=groups.setdefault(root(active[0]),(set(),[]))
            ids.update(active);rows.append(record)
        for ids,records in groups.values():
            ids=tuple(sorted(ids));indices={old:new for new,old in enumerate(ids)}
            local=FiniteConstraints()
            local.domains=[set(self.domains[i]) for i in ids]
            def translate(value):
                domain=self.domain(value)
                return next(iter(domain)) if len(domain)==1 else Choice(indices[value.index])
            for kind,values,data in records:
                if kind=='lookup':data={p:translate(v) for p,v in data.items()}
                local.constraints.append((kind,tuple(map(translate,values)),data))
            def token(value):return ('choice',value.index) if isinstance(value,Choice) else ('fixed',value)
            key=(tuple(tuple(sorted(d)) for d in local.domains),tuple(
                (kind,tuple(map(token,values)),tuple((p,token(v)) for p,v in sorted(data.items()))
                 if kind=='lookup' else data) for kind,values,data in local.constraints))
            yield ids,local,key

    def solve(self, timeout_ms):
        """Decide the reference-restricted components without enumerating mappings."""
        import z3
        deadline=time.perf_counter()+timeout_ms/1000
        self.propagate()
        stats=dict(components=0,cache_hits=0,solver_calls=0,
                   fixed_choices=sum(len(d)==1 for d in self.domains))
        if self.failed:return z3.unsat,None,stats
        assignments={i:min(d) for i,d in enumerate(self.domains)}
        for ids,component,key in self.components():
            stats['components']+=1
            if key in _component_cache:
                stats['cache_hits']+=1
                chosen=_component_cache[key];_component_cache.move_to_end(key)
                if chosen is None:return z3.unsat,None,stats
            else:
                solver=z3.SolverFor('QF_FD');component.emit(solver)
                remaining=deadline-time.perf_counter()
                if remaining<=0:return z3.unknown,None,stats
                solver.set(timeout=max(1,int(remaining*1000)))
                stats['solver_calls']+=1
                status=solver.check()
                if status==z3.unknown:return status,None,stats
                model=solver.model() if status==z3.sat else None
                def assignment(i):
                    value=component.expression(Choice(i))
                    return value if isinstance(value,int) else model.eval(value,model_completion=True).as_long()
                chosen=tuple(assignment(i) for i in range(len(ids))) if model is not None else None
                if len(repr(key))<=65536:
                    _component_cache[key]=chosen
                    if len(_component_cache)>256:_component_cache.popitem(last=False)
                if chosen is None:return z3.unsat,None,stats
            assignments.update(zip(ids,chosen))
        model=z3.Model()
        for i,domain in enumerate(self.domains):
            for value in domain:
                model.update_value(z3.Bool(f'q{i}_{value}'),z3.BoolVal(value==assignments[i]))
        return z3.sat,model,stats


@dataclass
class CompiledRelation:
    """Temporary symbols and factor provenance for one stored core relation."""
    system: FiniteConstraints
    relation: object
    images: dict
    discovery: dict
    source_action: tuple
    target_action: tuple
    program: tuple

    def _number(self, value, model):
        expression = self.system.expression(value[0])
        return expression if isinstance(expression, int) else model.eval(
            expression, model_completion=True).as_long()

    def choices(self, model):
        """Recover one discovery and both original-frame self permutations."""
        return dict(discovery={r: self._number(value, model)
                               for r, value in self.discovery.items()},
                    source_permutation=tuple(self._number(value, model) for value in self.source_action),
                    target_permutation=tuple(self._number(value, model) for value in self.target_action))

    def realize(self, model):
        mapping = self.relation.realize_permutations(**self.choices(model))
        images = {r: self._number(value, model) for r, value in self.images.items()}
        if mapping != images:
            raise ValueError('solver model disagrees with its relaxed relation choices')
        return mapping


def _compile_self(system, relation, label):
    """One correlated binary action, preserving elements and fixed anchors."""
    start = len(system.program)
    action = tuple(system.act([system.constant(i) for i in range(len(relation.atoms))],
                              relation.generators, label))
    return action, tuple(system.program[start:])


def _inverse_action(system, program, degree):
    """Invert the same selected factors, never introduce independent choices."""
    values = [system.constant(i) for i in range(degree)]
    for _, _, choice, factor in reversed(program):
        inverses = []
        for permutation in factor:
            inverse = [0] * len(permutation)
            for atom, image in enumerate(permutation):
                inverse[image] = atom
            inverses.append(tuple(inverse))
        values = system.apply_factor(values, tuple(inverses), choice)
    return values


def compile_discovery(system, relation, check=None):
    """Shared domains, occupancy, anchors and accepted weighted growth edges."""
    discovery = {r: system.symbol(domain) for r, domain in relation.domains}
    system.distinct(value[0] for value in discovery.values())
    for r, p in relation.anchors:
        system.restrict(discovery[r][0], {p})
    for edge in relation.accepted_edges:
        if check is not None:check()
        system.relation(tuple(discovery[r][0] for r in edge.atoms),
                        relation.discovery_pairs(edge))
    return discovery


def compile_relaxed_relation(system, relation, label, *, images_only=False):
    """Compile beta o D o alpha^-1 without expanding discovery maps or groups.

    Returned images use original endpoint indices. Both self actions use their
    compact structural factors, element labels and anchors.
    Only local edge tables are materialized; injectivity remains one shared
    constraint over the discovery choices. Image-only consumers may absorb
    actions already closed within discovery; callers querying action provenance
    retain the original factor variables.
    """
    discovery = compile_discovery(system, relation)
    source_absorbed, target_absorbed = relation.absorbed_self_actions if images_only else (False, False)
    source_action, source_program = (
        (tuple(system.constant(i) for i in range(len(relation.source_atoms))), ())
        if source_absorbed else _compile_self(system, relation.source, (label, 'source')))
    target_action, target_program = (
        (tuple(system.constant(i) for i in range(len(relation.target_atoms))), ())
        if target_absorbed else _compile_self(system, relation.target, (label, 'target')))
    inverse = _inverse_action(system, source_program, len(relation.source_atoms))
    original_targets = {i: system.constant(atom) for i, atom in enumerate(relation.target_atoms)}
    target_images = {atom: system.lookup(value, original_targets)
                     for atom, value in zip(relation.target_atoms, target_action)}
    discovery_table = {i: discovery[atom] for i, atom in enumerate(relation.source_atoms)}
    images = {atom: system.lookup(system.lookup(value, discovery_table), target_images)
              for atom, value in zip(relation.source_atoms, inverse)}
    return CompiledRelation(system, relation, images, discovery, source_action,
                            target_action, source_program + target_program)
