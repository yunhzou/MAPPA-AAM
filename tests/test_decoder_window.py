"""The window uses local constraints and only the selected sweep coupling."""
import pytest
import z3
from dataclasses import replace

from mappa import AAMProblem, decode_events
from mappa.mappa_core.post_processing.decoder.fragments import FragmentDecoder
from mappa.mappa_core.post_processing.decoder.records import EventDecodeConfig
from mappa.mappa_core.post_processing.final_branches import FinalBranchCatalogue
from mappa.mappa_core.utils.io.records import decoded_events_record
from growth_cases import saved_fragment
from test_event_interpretation import endpoint, saved, swap


def test_boundary_bound_rejects_composition_before_local_preparation(monkeypatch):
    ep = endpoint('CCCC', [(0, 1, 1.), (1, 2, 1.), (2, 3, 1.)])
    source = saved(ep, ep, dict(enumerate(range(4))))
    source.coupled_branches[0] = {tuple(((a,), (a,)) for a in range(4))}
    prepare = FragmentDecoder.fragment
    calls = []
    def checked(self, piece, *args):
        calls.append(piece.source_atoms)
        assert len(piece.source_atoms) == 4, 'fixed boundary already excludes singleton partition'
        return prepare(self, piece, *args)
    monkeypatch.setattr(FragmentDecoder, 'fragment', checked)
    result = decode_events(source, event_slack=0)
    record = decoded_events_record(result)
    assert result.minimum_events == result.event_window == 0
    assert len(record['tables']) == 1
    assert record['tables'][0]['fragments'] == ((0, 1, 2, 3),)
    assert len(result.local.fragments) == 1 and calls
    assert record['catalogue'] == source.to_record()


def test_sweep_resolves_only_its_endpoint_fragments_and_preserves_count_gaps(monkeypatch):
    monkeypatch.setattr(z3.Solver, 'check',
                        lambda *_: pytest.fail('Window construction must reuse prepared fragment counts'))
    r = endpoint('CCC', [(0, 1, 1.)])
    p = endpoint('CCC', [(0, 2, 1.)])
    fragment = saved_fragment(AAMProblem(r, p), source_atoms=(1, 2))
    source = saved(r, p, {0: 0, 1: 1, 2: 2}, fragments=(fragment,),
                   pairs=(((0,), (0,)), ((1, 2), (1, 2))))
    source.families[0] = replace(source.families[0], provenance=[dict(cuts=((0, 1),))])
    result = decode_events(source, event_slack=1)
    assert result.minimum_events == 0 and result.event_window == 1
    original = result.contributions(0)
    moved = result.contributions(0, [swap(fragment)])
    assert (original['total'], moved['total']) == (2, 0)
    assert not original['within_window'] and moved['within_window']
    assert [row.sum() for row in result.counts(0, [None, [swap(fragment)]])] == [2, 0]
    record = decoded_events_record(result)
    assert record['tables'][0]['event_window'] == 1
    assert record['minimum_scope'] == 'independent_fragment_model'
    assert record['event_counts_approximate'] and not record['global_minimum_proven']


def test_sweep_endpoint_choice_preserves_its_fragment_internal_cost(monkeypatch):
    from mappa.mappa_core.post_processing.decoder.local_counts import FragmentEventCounts
    monkeypatch.setattr(z3.Solver, 'check',
                        lambda *_: pytest.fail('Sweep assembly must use endpoint-conditioned counts'))
    r = endpoint('CCCC', [(0, 1, 1.), (1, 3, 1.)])
    p = endpoint('CCCC', [(0, 2, 1.), (1, 3, 1.)])
    fragment = saved_fragment(AAMProblem(r, p), source_atoms=(1, 2, 3), anchors=((3, 3),))
    source = saved(r, p, dict(enumerate(range(4))), fragments=(fragment,),
                   pairs=(((0,), (0,)), ((1, 2, 3), (1, 2, 3))))
    source.families[0] = replace(source.families[0], provenance=[dict(cuts=((0, 1),))])
    result = decode_events(source, event_slack=0)
    response = result.tables(0).fragment_responses[1][0]
    assert response.event_counts.minimum == 0
    bounds = response.event_counts.endpoint_bounds(1)
    assert bounds[1] == 0 and bounds[2] == 2
    # The same swap exchanges internal and sweep costs: (0 + 2) or (2 + 0).
    # Combining their separate zero minima would invent a nonexistent result.
    assert result.minimum_events == result.event_window == 2
    original = result.contributions(0)
    moved = result.contributions(0, [dict(discovery={1: 2, 2: 1, 3: 3})])
    assert original['total'] == moved['total'] == 2
    assert not original['internal'] and not moved['boundary']


def test_sweep_reuses_constant_fragment_costs_after_placement(monkeypatch):
    from mappa.mappa_core.core.symmetry.minimum import _Search

    # The 12-atom ring has many internal pair responses, but every rotation
    # costs zero. Only its attachment to the outside O can change the count.
    n = 12
    ring = [(i, 1 + i % n, 1.) for i in range(1, n + 1)]
    r = endpoint('O' + 'C' * n, ring + [(0, 1, 1.)])
    p = endpoint('O' + 'C' * n, ring + [(0, 2, 1.)])
    fragment = saved_fragment(AAMProblem(r, p), source_atoms=range(1, n + 1))
    source = saved(r, p, dict(enumerate(range(n + 1))), fragments=(fragment,),
        pairs=(((0,), (0,)), (tuple(range(1, n + 1)), tuple(range(1, n + 1)))))
    reflection = tuple(0 if i == 0 else 1 + (3 - i) % n for i in range(n + 1))
    source.families[0] = replace(source.families[0],
        actions=source.families[0].actions + (('group', (reflection,)),),
        provenance=[dict(cuts=((0, 1),))])
    solve = _Search.run
    def checked(self):
        if hasattr(self.data, 'actions'):
            pytest.fail('the sweep must project its two endpoints, not search complete choices')
        return solve(self)
    monkeypatch.setattr(_Search, 'run', checked)
    result = decode_events(source, event_slack=0)
    assert result.minimum_events == 0 and result.minimum_proven
    table = result.tables(0)
    local = table.fragment_responses[1][0].event_counts
    assert local.constant == 0 and not local.program.pairs
    assert len(table.program.pairs) == 2
    answer = result.select(0, 0).find()
    assert answer.status == 'found'
    assert result.contributions(0, answer.choices)['total'] == 0


def test_attained_contextual_bound_does_not_repeat_minimum_search(monkeypatch):
    from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents

    r = endpoint('CCC', [(0, 1, 1.), (1, 2, 1.)])
    p = endpoint('CCC', [(0, 1, 1.)])
    fragment = saved_fragment(AAMProblem(r, p))
    source = saved(r, p, dict(enumerate(range(3))), fragments=(fragment,),
        pairs=(((0, 1), (0, 1)), ((2,), (2,))))
    result = decode_events(source)
    counts = result.tables(0).fragment_responses[0][0].event_counts
    assert not counts.direct and counts.bounds == (0, 0)
    monkeypatch.setattr(OrderedSymmetryEvents, 'solve',
        lambda *a, **k: pytest.fail('a valid witness already attains the lower bound'))
    assert counts.minimum_result.complete and counts.minimum == 0


def test_weight_preserving_placement_reuses_a_variable_fragment_minimum(monkeypatch):
    from mappa.mappa_core.core.symmetry.program import OrderedSymmetryEvents

    source_weights, target_weights = (2.,2.,1.,1.,1.,1.), (2.,1.,2.,1.,1.,1.)
    def rings(values):
        return endpoint('C'*12, [(offset+i, offset+(i+1)%6, weight)
            for offset in (0,6) for i,weight in enumerate(values)])
    r, p = rings(source_weights), rings(target_weights)
    fragment = saved_fragment(AAMProblem(r,p), source_atoms=range(6),
                              accepted_edges=tuple((i,(i+1)%6) for i in range(6)))
    source = saved(r,p,dict(enumerate(range(6))),fragments=(fragment,),
                   pairs=((tuple(range(6)),tuple(range(6))),))
    source.families[0] = replace(source.families[0], actions=source.families[0].actions +
                                (('group', (tuple((i+6)%12 for i in range(12)),)),))
    monkeypatch.setattr(OrderedSymmetryEvents, 'solve',
        lambda *a, **k: pytest.fail('cost-preserving placements repeated the prepared fragment minimum'))
    decoded = decode_events(source)
    counts = decoded.tables(0).fragment_responses[0][0].event_counts
    assert counts.program.pairs and not counts.direct
    assert counts.minimum == 2 and counts.minimum_proven
    assert decoded.contributions(0, counts.minimum_choices)['total'] == 2


def test_minimum_at_initial_upper_endpoint_is_retained():
    r, p = endpoint('CCC', [(0, 1, 1.)]), endpoint('CCC', [(0, 2, 1.)])
    source = saved(r, p, {0: 0, 1: 1, 2: 2},
                   pairs=(((0,), (0,)), ((1,), (1,)), ((2,), (2,))))
    source.families[0] = replace(source.families[0], provenance=[dict(cuts=((0, 1),))])
    result = decode_events(source, event_slack=0)
    assert result.tables(0).boundary_events.lower == 0
    assert result.minimum_events == result.event_window == 2
    assert result.window.contexts == ((0, 0, 0),)
    assert result.contributions(0)['within_window']


def test_initial_sweep_window_uses_local_minima_before_the_reporting_window(monkeypatch):
    from mappa.mappa_core.post_processing.decoder.window import EventWindow

    # The sweep can contribute 0 or 2. A separate fixed path contributes
    # 4, 0, 2 or 6 events, although every composition has boundary lower bound 0.
    r = endpoint('CCCCCCC', [(0, 1, 1.), (3, 4, 1.), (4, 5, 1.), (5, 6, 1.)])
    p = endpoint('CCCCCCC', [(0, 2, 1.), (3, 4, 1.), (4, 5, 1.), (5, 6, 1.)])
    problem = AAMProblem(r, p)
    movable = saved_fragment(problem, source_atoms=(1, 2))
    source = None
    for images in ((3, 5, 4, 6), (3, 4, 5, 6), (3, 4, 6, 5), (4, 6, 3, 5)):
        path_mapping = tuple(zip(range(3, 7), images))
        fixed = saved_fragment(problem, source_atoms=range(3, 7),
                               witness=path_mapping, anchors=path_mapping)
        mapping = {0: 0, 1: 1, 2: 2, **dict(path_mapping)}
        row = saved(r, p, mapping, fragments=(movable, fixed),
                    pairs=(((0,), (0,)), ((1, 2), (1, 2)), ((3, 4, 5, 6), (3, 4, 5, 6))))
        row.families[0] = replace(row.families[0], provenance=[dict(cuts=((0, 1),))])
        if source is None:
            source = row
        else:
            source.families.extend(row.families)
            next(iter(source.branches.values())).add(len(source.families) - 1)

    result = decode_events(source, event_slack=4)
    assert result.minimum_events == 0 and result.event_window == 4
    assert result.window.initial_lower == 0
    assert result.window.initial_maximum <= 2
    assert all(result.window.bounds[key] <= 2 for key in result.window.results)
    assert result.window.contexts == ((0, 0, 0), (1, 0, 0), (2, 0, 0))
    # The lower-4 composition was excluded from minimum finding, but belongs
    # to the larger reporting window. The lower-6 composition belongs to neither.
    assert [result.contributions(fid, [swap(movable), None])['total']
            for fid in range(4)] == [4, 0, 2, 6]
    assert [result.contributions(fid)['total'] for fid in range(4)] == [6, 2, 4, 8]
    shared = [response for response in result.local.fragments.values()
              if response.piece.source_atoms == (1, 2)]
    assert len(shared) == 1
    assert shared[0].event_counts.minimum == 0


def test_incompatible_local_minima_do_not_invent_a_feasible_initial_window():
    # This manually partitioned whole-fragment relation lets both pieces
    # minimize onto the same target CH3. Distinct sweep endpoints prevent it.
    r = endpoint('CHHHCHHH', [(0, 1, 1.), (0, 2, 1.), (0, 3, 1.),
                             (4, 5, 1.), (4, 6, 1.), (4, 7, 1.), (0, 4, 1.)])
    p = endpoint('CHHHCHHH', [(0, 1, 1.), (0, 2, 1.), (0, 3, 1.)])
    fragment = saved_fragment(AAMProblem(r, p))
    source = saved(r, p, dict(enumerate(range(8))), fragments=(fragment,),
                   pairs=(((0, 1, 2, 3), (0, 1, 2, 3)), ((4, 5, 6, 7), (4, 5, 6, 7))))
    source.families[0] = replace(source.families[0], provenance=[dict(cuts=((0, 4),))])
    result = decode_events(source)
    table = result.tables(0)
    assert table.boundary_events.lower == 0
    assert [response.event_counts.minimum for response, _ in table.fragment_responses] == [0, 0]
    assert result.minimum_events == 4
    answer = result.minimum(0)
    assert answer.complete and answer.event_count == 4
    assert len(set(dict(answer.witness).values())) == 8
    assert result.contributions(0, answer.choices)['total'] == 4


def test_reporting_window_reuses_prepared_fragment_counts_without_solver_queries(monkeypatch):
    ep = endpoint('CHHHON', [(0, 1, 1.), (0, 2, 2.), (0, 3, 3.), (0, 4, 1.), (0, 5, 1.)])
    source = saved(ep, ep, dict(enumerate(range(6))))
    mapping = {0: 0, 1: 2, 2: 1, 3: 3}
    fragment = saved_fragment(AAMProblem(ep, ep), source_atoms=range(4), witness=mapping)
    partial = saved(ep, ep, mapping, fragments=(fragment,))
    source.families.extend(partial.families)
    source.branches.update({key: {1} for key in partial.branches})

    monkeypatch.setattr(z3.Solver, 'check',
                        lambda *_: pytest.fail('Reporting must reuse prepared fragment responses'))
    result = decode_events(source, event_slack=6)
    record = decoded_events_record(result)
    assert result.minimum_events == 0 and result.event_window == 6
    # Four ordinary boundary events exclude the partial star from L+2.
    # Its local count remains conditional between 0 and 2 in the M+6 output.
    assert result.tables(1).count_bounds == (4, 6)
    assert result.window.bounds[(1, 0, 0)] == 4
    assert len(record['tables']) == 2
    assert record['catalogue'] == source.to_record()


def test_separate_sweep_histories_do_not_become_a_two_cut_query():
    ep = endpoint('CCC', [(0, 1, 1.), (1, 2, 1.)])
    source = saved(ep, ep, {0: 0, 1: 1, 2: 2},
                   pairs=(((0,), (0,)), ((1,), (1,)), ((2,), (2,))))
    source.families[0] = replace(source.families[0], provenance=[dict(cuts=((0, 1),)), dict(cuts=((1, 2),))])
    result = decode_events(source)
    assert result.cut_contexts == (((((0, 1),), ((1, 2),)),),)
    assert result.minimum_events == 2
    assert result.tables(0, cut=0).cut_bonds == frozenset(((0, 1),))
    assert result.tables(0, cut=1).cut_bonds == frozenset(((1, 2),))
    assert len(decoded_events_record(result)['tables']) == 2


def test_ordinary_boundary_approximation_is_separate_from_physical_report():
    ep = endpoint('CC', [(0, 1, 1.)])
    source = saved(ep, ep, {0: 0, 1: 1}, pairs=(((0,), (0,)), ((1,), (1,))))
    result = decode_events(source)
    assert result.minimum_events == result.contributions(0)['total'] == 2
    assert result.evaluate(0)['total'] == 0
    assert result.realize(0) == {0: 0, 1: 1}


def test_empty_catalogue_has_no_minimum_or_window():
    ep = endpoint('C', [])
    result = decode_events(FinalBranchCatalogue(AAMProblem(ep, ep)))
    assert result.minimum_events is None and result.event_window is None
    assert decoded_events_record(result)['tables'] == []


def test_fixed_sweep_response_needs_no_joint_fragment_query(monkeypatch):
    r, p = endpoint('CCCC', [(0, 2, 1.)]), endpoint('CCCC', [])
    problem = AAMProblem(r, p)
    fragments = tuple(saved_fragment(problem, source_atoms=atoms) for atoms in ((0, 1), (2, 3)))
    source = saved(r, p, dict(enumerate(range(4))), fragments=fragments,
                   pairs=(((0, 1), (0, 1)), ((2, 3), (2, 3))))
    source.families[0] = replace(source.families[0], provenance=[dict(cuts=((0, 2),))])
    def forbidden(*args):
        pytest.fail('the sweep costs one event for every allowed endpoint pair')
    monkeypatch.setattr(z3.Solver, 'check', forbidden)
    result = decode_events(source)
    assert result.minimum_events == 1
    for choice in ([None, None], [swap(fragment) for fragment in fragments]):
        assert result.contributions(0, choice)['total'] == 1


def test_partition_provenance_keeps_its_own_sweep_context():
    ep = endpoint('CCC', [(0, 1, 1.), (1, 2, 1.)])
    source = saved(ep, ep, {0: 0, 1: 1, 2: 2},
                   pairs=(((0,), (0,)), ((1, 2), (1, 2))))
    source.coupled_branches[0] = {(((0, 1), (0, 1)), ((2,), (2,)))}
    source.families[0] = replace(source.families[0], provenance=[
        dict(cuts=((0, 1),), fragment_partition=((0,), (1, 2))),
        dict(cuts=((1, 2),), fragment_partition=((0, 1), (2,)))])
    result = decode_events(source)
    assert result.cut_contexts == (((((0, 1),),), (((1, 2),),)),)
    assert len(decoded_events_record(result)['tables']) == 2


@pytest.mark.parametrize('slack', [-1, 1.5, None, True])
def test_window_has_one_nonnegative_integer_setting(slack):
    with pytest.raises(ValueError, match='event_slack'):
        EventDecodeConfig(event_slack=slack)
