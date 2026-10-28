from dataclasses import replace
import json
import pickle
import pytest
from mappa.execution import Task, TaskQueue
from rebenchmark.decoder_inventory import Inventory, execute


def task(seed=0):
    return Task(f'00000/cut/00000/seed/{seed:02d}', 'rebenchmark.decoding:decode_seed',
                dict(task_id=0, cut=0, seed=seed, slack=5), case='00000')


def inventory(tmp_path):
    return Inventory(tmp_path/'inventory.sqlite', [dict(task=0, cuts=1, expected=2)], seeds=2)


def test_missing_and_execution_completion_do_not_prove_windows(tmp_path):
    i = inventory(tmp_path)
    t = task()
    i.replace(t.key, 'complete', {'blob': 'output'}, 'original', {'missing'}, t)
    assert i.summary()['states'] == {'complete': 1, 'missing': 1}
    assert not i.summary()['scientific_audit_complete']


def test_continuation_preserves_science_and_allows_remaining_budget(tmp_path):
    i = inventory(tmp_path); t = task()
    i.replace(t.key, 'pending', None, 'original', {'missing'}, t)
    changed = replace(t, payload=dict(t.payload, slack=6))
    with pytest.raises(ValueError, match='Scientific'):
        i.replace(t.key, 'complete', {'blob': 'wrong'}, 'continuation', {'pending'}, changed)
    i.replace(t.key, 'complete', {'blob': 'saved'}, 'continuation', {'pending'}, replace(t, case_seconds=50))
    with pytest.raises(ValueError, match='duplicate'):
        i.replace(t.key, 'complete', {'blob': 'again'}, 'repair', {'error'}, t)


def test_recovery_and_pilots_are_disjoint(tmp_path):
    i = inventory(tmp_path); t = task()
    i.replace(t.key, 'lost', None, 'original', {'missing'}, t)
    i.recovered({t.key: dict(key=t.key, blob='saved')}, 'recovered')
    p = task(1)
    i.pilots({p.key: dict(state='complete', resource={'blob': 'pilot'})}, 'pilots')
    assert i.summary()['states'] == {'complete': 2}
    with pytest.raises(ValueError, match='duplicate'):
        i.recovered({t.key: dict(key=t.key, blob='duplicate')}, 'second_recovery')


def test_wrong_task_identity_and_missing_output_rejected(tmp_path):
    i = inventory(tmp_path); t = task()
    with pytest.raises(ValueError, match='identity'):
        i.replace(t.key, 'pending', None, 'original', {'missing'}, replace(t, key='wrong'))
    with pytest.raises(ValueError, match='lacks'):
        i.replace(t.key, 'complete', {}, 'original', {'missing'}, t)


def test_seed_grid_must_be_complete(tmp_path):
    with pytest.raises(ValueError, match='grid'):
        Inventory(tmp_path/'bad.sqlite', [dict(task=0,cuts=1,expected=1)], seeds=2)


def test_full_overlay_accounting(tmp_path):
    def write_json(name, value):
        p = tmp_path/name
        p.write_text(json.dumps(value))
        return str(p)

    def queue(name, rows):
        p = tmp_path/name
        q = TaskQueue(p, name)
        q.add([t for t, _ in rows])
        with q.db:
            for t, state in rows:
                resource = {'blob': name+t.key} if state == 'complete' else {'error': 'diagnosed'}
                q.db.execute('UPDATE tasks SET state=?,result=? WHERE key=?', (state,pickle.dumps(resource),t.key))
        q.close()
        return str(p)

    original = queue('original.sqlite', [(task(s), state) for s,state in enumerate(['complete','pending','lost','error'])])
    continuation = queue('continuation.sqlite', [(replace(task(1),case_seconds=50), 'complete')])
    repair = queue('repair.sqlite', [(replace(task(3),case_seconds=20), 'complete')])
    plan = dict(seeds=5, original=original, continuation=continuation,
        expected=write_json('expected.json', dict(total=5,directions=[dict(task=0,cuts=1,expected=5)])),
        pilots=write_json('pilots.json', {task(4).key: dict(state='complete',resource={'blob':'pilot'})}),
        recovery=write_json('recovery.json', dict(missing=[],duplicates=[],recovered={task(2).key: dict(key=task(2).key,blob='recovered')})),
        repairs=[dict(name='repair',queue=repair,selection=write_json('selection.json', dict(selected=[dict(key=task(3).key)])))])
    summary = execute(plan, tmp_path/'inventory.sqlite')
    assert summary['states'] == {'complete': 5}
    assert summary['sources'] == dict(original=4,pilots=1,startup_recovery=1,continuation=1,repair=1)
    assert summary['scientific_audit_complete'] is False
