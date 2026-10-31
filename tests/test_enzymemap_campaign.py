"""Fresh benchmark orchestration uses public stages and bounded archived history."""
from dataclasses import asdict
import json
from pathlib import Path
import sqlite3
import sys

import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from enzymemap import run
from mappa import AAMSearchConfig
from mappa.mappa_core.core.scheduler import Task
from mappa.mappa_core.core.scheduler.storage import BundleStore
from mappa.mappa_core.utils.io.records import endpoint_record
from test_reused_aam_api import problem


def prepared_index(tmp_path, *, corrupt=False):
    store=BundleStore(tmp_path/'prepared')
    case=problem()
    reference=store.put(dict(mapping=[(i,i) for i in range(4)]))
    rows=[]
    for reverse in (False,True):
        left,right=(case.product,case.reactant) if reverse else (case.reactant,case.product)
        source=store.put(dict(name="test",reactant=endpoint_record(left),product=endpoint_record(right)))
        saved=store.put([Task('retired','retired.must_never_execute:run',dict(input=source))])
        row=dict(case=0,reversed=reverse,blob=asdict(saved),reference=asdict(reference))
        if corrupt: row['blob']['sha256']='0'*64
        rows.append(row)
    path=tmp_path/'index.json'
    path.write_text(json.dumps(dict(config=asdict(AAMSearchConfig(seed_count=1,branch_limit=10)),directions=rows)))
    return path


def test_full_passes_and_resume_preserve_results_without_old_operations(tmp_path):
    index=prepared_index(tmp_path)
    output=tmp_path/'result'
    status=run(index,output,local=True,workers=2,window=2)
    assert status['finished']
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute("SELECT count(*) FROM results WHERE state='complete'").fetchone()[0] == 8
        record=db.execute("SELECT reference FROM results WHERE stage='verification'").fetchone()[0]
    from mappa.mappa_core.core.scheduler.bundles import Blob
    assert BundleStore(output/'bundles').load(Blob(**json.loads(record)))['status']=='recovered'
    with sqlite3.connect(output/'bundles/queue.sqlite') as db:
        assert db.execute('SELECT count(*) FROM tasks').fetchone()[0] == 0
    before=len(list((output/'bundles').glob('*.bundle')))
    assert run(index,output,local=True,workers=1,resume=True)['finished']
    assert len(list((output/'bundles').glob('*.bundle'))) == before


def test_input_failures_remain_in_denominator_and_stop_next_pass(tmp_path):
    index=prepared_index(tmp_path,corrupt=True)
    output=tmp_path/'result'
    with pytest.raises(RuntimeError,match='matching contains failed'):
        run(index,output,local=True,workers=1)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT stage,state,count(*) FROM results GROUP BY stage,state').fetchall()==[('matching','error',2)]


def test_retry_ledger_retains_original_failure_and_audit(tmp_path):
    from enzymemap import Ledger
    output=tmp_path/'result'; output.mkdir()
    scratch=tmp_path/'scratch'; scratch.mkdir()
    store=BundleStore(output/'bundles')
    audit=store.put({'attempts': ['original memory failure']})
    ledger=Ledger(output,scratch,resume=False)
    ledger.put('case','matching','error',audit=audit,error='memory_limit')
    assert ledger.retry_audit('case','matching')==audit
    ledger.retire_failure('case','matching','diagnosed native export defect')
    ledger.put('case','matching','complete',reference=store.put('replacement'))
    ledger.checkpoint(); ledger.db.close()
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT state FROM results').fetchone()==('complete',)
        reason,raw,error=db.execute('SELECT reason,audit,error FROM retries').fetchone()
        assert reason=='diagnosed native export defect'
        assert json.loads(raw)==asdict(audit)
        assert json.loads(error)=='memory_limit'


def test_resource_retry_recovers_old_timeout_but_keeps_same_memory_failure(tmp_path):
    from mappa.mappa_core.core.scheduler.bundles import Blob
    index = prepared_index(tmp_path)
    output = tmp_path/'result'
    run(index, output, local=True, workers=2, stop_after='matching')
    store = BundleStore(output/'bundles')
    originals = {}
    with sqlite3.connect(output/'ledger.sqlite') as db:
        rows = db.execute('SELECT case_name,audit FROM results ORDER BY case_name').fetchall()
        for (case, raw), memory in zip(rows, (2*2**30, 4*2**30)):
            audit = store.load(Blob(**json.loads(raw)))
            originals[case] = audit
            task = list(audit['tasks'][-1])
            task[3], task[4] = memory, 'timeout' if memory == 2*2**30 else 'memory_limit'
            audit['tasks'][-1] = tuple(task)
            budget = list(audit['case']); budget[3] = budget[1]
            audit['case'] = tuple(budget)
            db.execute("UPDATE results SET state='error',reference=NULL,audit=? WHERE case_name=?",
                (json.dumps(asdict(store.put(audit))), case))
    with pytest.raises(RuntimeError, match='matching contains failed'):
        run(index, output, local=True, workers=2, resume=True, stop_after='matching',
            retry_memory_limits=True, case_seconds=60.)
    with sqlite3.connect(output/'ledger.sqlite') as db:
        assert db.execute('SELECT case_name,state FROM results ORDER BY case_name').fetchall() == [
            ('00000/forward', 'complete'), ('00000/reverse', 'error')]
        assert db.execute('SELECT case_name FROM retries').fetchall() == [('00000/forward',)]
        raw, = db.execute("SELECT audit FROM results WHERE case_name='00000/forward'").fetchone()
    recovered = store.load(Blob(**json.loads(raw)))
    original = originals['00000/forward']
    assert len(recovered['attempts']) == len(original['attempts'])+1
    assert recovered['tasks'][:-1] == original['tasks'][:-1]
    assert recovered['retries'][0][4:] == (60., 4*2**30)


def test_batches_are_disjoint_reuse_results_and_stop_before_decoding(tmp_path):
    from prepare_batches import prepare
    index=prepared_index(tmp_path)
    data=json.loads(index.read_text())
    data['directions'] += [dict(row,case=1) for row in data['directions']]
    index.write_text(json.dumps(data))
    baseline=tmp_path/'baseline'
    assert run(index,baseline,local=True,workers=2,stop_after='matching')['pass_finished']=='matching'
    prepare(index,baseline,tmp_path/'batches',count=2,
            additional_ledgers=[baseline/'ledger.sqlite'])
    for batch in range(2):
        root=tmp_path/'batches'/str(batch)
        assert {r['case'] for r in json.loads((root/'index.json').read_text())['directions']}=={batch}
        with sqlite3.connect(root/'result/ledger.sqlite') as db:
            assert db.execute('SELECT stage,count(*) FROM results GROUP BY stage').fetchall()==[('matching',2)]
        alias=tmp_path/f'alias-{batch}'
        alias.symlink_to(root, target_is_directory=True)
        assert run(alias/'index.json',alias/'result',local=True,workers=1,resume=True)['finished']
