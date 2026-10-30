"""Immutable dataset preparation and compatibility checks for reusable inputs."""
from pathlib import Path
import tempfile
from mappa.mappa_core.core.scheduler.search import plan_search
from .operations import conversion_tasks, reference_tasks
from .evidence import read, save, sha


def _initialize_planner(output, scratch):
    from mappa.mappa_core.core.scheduler import BundleWriter
    global _planning_writer
    _planning_writer = BundleWriter(Path(output)/'inputs', scratch)


def _prepare_record(payload):
    """Plan one reaction, publishing only bounded bundles from local scratch."""
    import hashlib
    import json
    import numpy as np
    from collections import Counter
    from dataclasses import asdict
    from mappa import AAMSearchConfig, AAMProblem
    from mappa.mappa_core.utils.io.records import endpoint_record, problem_from_record
    plan, item, output = payload
    case = item['index']
    if plan['dataset'] == 'EnzymeMap':
        from golden_evaluation import prepare
        problem, input_features, _ = prepare(item['rxn'])
        problem = AAMProblem(problem.reactant, problem.product, f'enzymemap_{case}')
        reference_problem, features, mapping = prepare(item['ground_truth'])
        same_index = all(a.elements == b.elements and np.array_equal(a.wbo, b.wbo)
            for a, b in zip((problem.reactant, problem.product),
                            (reference_problem.reactant, reference_problem.product)))
        raw = dict(input=dict(name=problem.name,
            reactant=endpoint_record(problem.reactant, array_lists=True),
            product=endpoint_record(problem.product, array_lists=True)),
            reference=dict(features=features, mapping=sorted(mapping.items()),
                mapped_reaction=item['ground_truth'], input_reaction=item['rxn'],
                input_features=input_features, input_reference_index_agrees=same_index,
                original_id=item['original_id']))
    else:
        folder = Path(plan['input_manifest']).parent/'inputs'/str(case)
        raw = {}
        for name in ('input', 'reference'):
            data = (folder/f'{name}.json').read_bytes()
            if hashlib.sha256(data).hexdigest() != item[f'{name}_sha256']:
                raise ValueError(f'Input {case} {name} hash differs')
            raw[name] = json.loads(data)
        problem = problem_from_record(raw['input'])
    config = AAMSearchConfig(**plan['config'])
    writer = _planning_writer
    reference = asdict(writer.put(raw['reference']))
    rows = []
    for reverse in plan['directions']:
        direction = 2*case+int(reverse)
        directed = AAMProblem(problem.product,problem.reactant,problem.name) if reverse else problem
        tasks, cuts = plan_search(f'{direction:05d}', directed, config, writer,
            seed_batch=plan['seed_batch'], case_seconds=plan['case_seconds'], memory_bytes=plan['memory_bytes'])
        tasks += list(conversion_tasks(f'{direction:05d}',tasks,cuts,slack=plan['slack']))
        if plan.get('reference_checks', False):
            tasks += reference_tasks(f'{direction:05d}', tasks, cuts, reference, reverse)
        order = {'search_batch':0,'finalize_cut':1,'verify_seed':2,'convert_seed':3}
        tasks.sort(key=lambda t: (t.key.split('/')[2], order[t.operation.rsplit(':',1)[1]], t.key))
        rows.append(dict(direction=direction,case=case,reversed=reverse,cuts=len(cuts),
            tasks=len(tasks),counts=dict(Counter(t.operation.rsplit(':',1)[1] for t in tasks)),
            blob=asdict(writer.put(tasks)),reference=reference,
            input_sha256=hashlib.sha256(json.dumps(raw['input'],sort_keys=True).encode()).hexdigest(),
            reference_sha256=hashlib.sha256(json.dumps(raw['reference'],sort_keys=True).encode()).hexdigest()))
    return rows


def prepare_campaign(out):
    """Freeze a dataset's explicit cut/seed grid, with parallel input planning."""
    import csv
    import gzip
    import hashlib
    import json
    from collections import Counter
    from concurrent.futures import ProcessPoolExecutor
    from dataclasses import asdict
    from mappa import AAMSearchConfig
    from .evidence import save
    plan = json.loads((out/'plan.json').read_text())
    manifest_path = Path(plan['input_manifest'])
    if hashlib.sha256(manifest_path.read_bytes()).hexdigest() != plan['input_manifest_sha256']:
        raise ValueError('Input manifest changed')
    manifest = json.loads(manifest_path.read_text())
    plan.setdefault('dataset', 'Golden')
    plan.setdefault('directions', [False, True])
    if plan['dataset'] == 'EnzymeMap':
        dataset = Path(plan['dataset_path'])
        if hashlib.sha256(dataset.read_bytes()).hexdigest() != manifest['dataset_sha256']:
            raise ValueError('EnzymeMap dataset changed')
        with gzip.open(dataset,'rt') as stream:
            items = [dict(row,index=int(row['r_id'])) for row in csv.DictReader(stream)]
        if len(items) != manifest['cases'] or [r['index'] for r in items] != list(range(manifest['cases'])):
            raise ValueError('EnzymeMap row identities differ')
    elif plan['dataset'] == 'Golden':
        items = manifest['inputs']
        if len(items) != 1851 or len(manifest['tasks']) != 3702:
            raise ValueError('Golden input grid is incomplete')
    else:
        raise ValueError('Unknown benchmark dataset')
    config = AAMSearchConfig(**plan['config'])
    if (out/'index.json').exists():
        raise FileExistsError('The frozen plan already exists')
    if 'selected_cases' in plan:
        selected = set(plan['selected_cases'])
        items = [r for r in items if r['index'] in selected]
        if len(items) != len(selected):
            raise ValueError('Selected cases are missing')
    rows, totals = [], Counter()
    import os
    with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as scratch:
        with ProcessPoolExecutor(max_workers=plan.get('preparation_workers',1),
                initializer=_initialize_planner,initargs=(str(out),scratch)) as pool:
            for done, prepared in enumerate(pool.map(_prepare_record,
                    ((plan,item,str(out)) for item in items), chunksize=1),1):
                rows.extend(prepared)
                for row in prepared:
                    totals.update(row['counts'])
                if done % 100 == 0:
                    save(out/'preparation.json',dict(cases=done,tasks=sum(totals.values())))
    save(out/'index.json',dict(schema='mappa.fresh_benchmark/v1',dataset=plan['dataset'],
        config=asdict(config),directions=rows,counts=dict(totals),tasks=sum(totals.values()),
        cases=len(items),input_manifest_sha256=plan['input_manifest_sha256'],
        preparation=specification(plan)))
    save(out/'preparation.json',dict(finished=True,cases=len(items),tasks=sum(totals.values())))


def load_plan(folder):
    from mappa import AAMProblem,MolecularEndpoint,AAMSearchConfig,AAMSearchPlan
    raw=read(Path(folder)/'input.json');task=read(Path(folder)/'task.json')
    original=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')),raw['name'])
    reverse=task['direction']=='P_to_R'
    directed=AAMProblem(original.product,original.reactant,original.name) if reverse else original
    return AAMSearchPlan(original,directed,AAMSearchConfig(**task['config']),reverse)


def specification(plan):
    """Settings that define prepared tasks; runtime placement is deliberately absent."""
    from dataclasses import asdict
    from mappa import AAMSearchConfig
    import json
    return json.loads(json.dumps(dict(dataset=plan['dataset'], input_manifest_sha256=plan['input_manifest_sha256'],
        config=asdict(AAMSearchConfig(**plan['config'])), directions=plan['directions'],
        selected_cases=plan.get('selected_cases'), seed_batch=plan['seed_batch'],
        case_seconds=plan['case_seconds'], memory_bytes=plan['memory_bytes'],
        slack=plan['slack'], reference_checks=plan['reference_checks'],
        task_schema='mappa.benchmark.operations/v1')))


def validate(plan):
    """An explicit artifact is reusable only with the same task semantics."""
    path = Path(plan['index'])
    if sha(path) != plan['index_sha256']:
        raise ValueError('Prepared input index changed')
    index = read(path)
    if index['preparation'] != specification(plan):
        raise ValueError('Prepared inputs differ from requested settings; prepare a new artifact')
    if index['cases'] != plan['expected_cases']:
        raise ValueError('Prepared dataset denominator differs')
    return index


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    prepare_campaign(parser.parse_args().output)
