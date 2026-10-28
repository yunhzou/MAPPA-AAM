"""Finite reference constraints, reduced before constructing the exact solver query.

Each permutation factor has one shared choice variable. Local propagation only
removes values without support; the solver retains all correlations afterward.
"""
from collections import deque
from dataclasses import dataclass


@dataclass(frozen=True)
class Choice:
    index: int


class ReferenceConstraints:
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
        from .family_query import group_factors
        for factor in group_factors(tuple(tuple(g) for g in generators)):
            active=set().union(*(v[1] for v in values))
            if not any(p<len(factor[0]) and any(g[p]!=p for g in factor) for p in active):continue
            choice=self.symbol(range(len(factor)))
            transformed=[]
            for value in values:
                rows=tuple((p,c,g[p] if p<len(g) else p) for p in value[1] for c,g in enumerate(factor))
                if all(p==q for p,c,q in rows):transformed.append(value);continue
                output=self.symbol(q for p,c,q in rows)
                self.relation((value[0],choice[0],output[0]),rows)
                transformed.append(output)
            values=transformed
            self.program.append(('group',label,choice[0],factor))
        return values

    def propagate(self):
        neighbors=[set() for _ in self.domains]
        for n,(kind,values,data) in enumerate(self.constraints):
            for value in (*values,*(data.values() if kind=='lookup' else ())):
                if isinstance(value,Choice):neighbors[value.index].add(n)
        queue=deque(range(len(self.constraints)));queued=set(queue)
        while queue and not self.failed:
            n=queue.popleft();queued.remove(n)
            kind,values,data=self.constraints[n];changed=set()
            if kind=='table':
                valid=[]
                for row in data:
                    assignments={}
                    for value,item in zip(values,row):
                        if item not in self.domain(value) or (value in assignments and assignments[value]!=item):break
                        assignments[value]=item
                    else:valid.append(row)
                if not valid:self.failed=True;break
                for i,value in enumerate(values):changed |= self.restrict(value,{r[i] for r in valid})
                self.constraints[n]=(kind,values,tuple(valid))
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
        import z3
        from itertools import product
        self.propagate()
        if self.failed:solver.add(False);return
        for i,domain in enumerate(self.domains):
            if len(domain)>1:
                solver.add(z3.PbEq([(self.predicate(Choice(i),v),1) for v in sorted(domain)],1))
        for kind,values,data in self.constraints:
            if kind=='table':
                # A local implication per input tuple avoids asking the solver
                # to select from a large disjunction of complete assignments.
                supported={}
                for row in data:supported.setdefault(row[:-1],set()).add(row[-1])
                for inputs in product(*(sorted(self.domain(v)) for v in values[:-1])):
                    premise=[z3.Not(self.predicate(v,x)) for v,x in zip(values[:-1],inputs)]
                    outputs=[self.predicate(values[-1],x) for x in sorted(supported.get(inputs,()))]
                    solver.add(z3.Or(*premise,*outputs))
            elif kind=='lookup':
                selector,output=values
                for item in sorted(self.domain(selector)):
                    entry=data[item]
                    for result in sorted(self.domain(output)|self.domain(entry)):
                        solver.add(z3.Implies(self.predicate(selector,item),
                            self.predicate(output,result)==self.predicate(entry,result)))
            else:
                for item in sorted(set().union(*(self.domain(v) for v in values))):
                    flags=[self.predicate(v,item) for v in values]
                    fixed=sum(flag is True for flag in flags)
                    choices=[flag for flag in flags if not isinstance(flag,bool)]
                    if fixed>1:solver.add(False)
                    elif choices:solver.add(z3.PbLe([(flag,1) for flag in choices],1-fixed))
