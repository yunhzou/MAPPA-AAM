from dataclasses import replace

from mappa import AAMSearchConfig, search_aam
from mappa.mappa_core.core.search.graph import AAMSearchGraph
from native_cases import chain_problem as problem


def old_combine(graphs):
    contexts, roots, states, edges, stops = [], [], [], [], []
    for graph in graphs:
        c,s,t = len(contexts),len(states),len(edges)
        contexts.extend(graph.contexts)
        roots.extend(r+s for r in graph.roots)
        states.extend(replace(v,id=v.id+s,context=v.context+c) for v in graph.states)
        edges.extend(replace(e,id=e.id+t,source=e.source+s,target=e.target+s) for e in graph.transitions)
        stops.extend(replace(v, state=None if v.state is None else v.state+s,
                             context=graph.stop_context(v)+c if v.state is None else None)
                     for v in graph.stops)
    return AAMSearchGraph(tuple(contexts),tuple(roots),tuple(states),tuple(edges),tuple(stops))


def test_direct_offset_construction_preserves_every_field():
    graphs = [search_aam(problem(n), AAMSearchConfig(seed_count=2,branch_limit=cap)).graph
              for n,cap in [(4,1),(6,100),(3,3)]]
    for values in ([],graphs[:1],graphs,graphs+graphs):
        assert AAMSearchGraph.combine(values) == old_combine(values)


def test_collection_references_prepared_fragments_without_reconstruction(monkeypatch):
    import pickle
    from mappa.mappa_core.core.matching.representation.records import Fragment
    graph = search_aam(problem(4), AAMSearchConfig(seed_count=1, sweep_cuts=False)).graph
    edge = next(edge for edge in graph.transitions if edge.match is not None)
    def reconstruct(*args, **kwargs):
        raise AssertionError('collection must reference completed fragments directly')
    monkeypatch.setattr(Fragment, 'from_record', reconstruct)
    combined = AAMSearchGraph.combine((graph,) * 100)
    references = [combined.transitions[i * len(graph.transitions) + edge.id].match
                  for i in range(100)]
    assert all(fragment is edge.match for fragment in references)
    restored = pickle.loads(pickle.dumps(combined))
    references = [restored.transitions[i * len(graph.transitions) + edge.id].match
                  for i in range(100)]
    assert all(fragment is references[0] for fragment in references)
    assert restored == combined
