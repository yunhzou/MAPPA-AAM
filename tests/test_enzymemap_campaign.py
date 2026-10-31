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
