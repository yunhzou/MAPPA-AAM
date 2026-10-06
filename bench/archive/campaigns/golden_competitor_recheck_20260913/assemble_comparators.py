"""Audit persisted competitor records and package the verified paper evidence."""
from pathlib import Path
import json,hashlib,gzip,shutil,ast
from collections import Counter
V=Path(__file__).resolve().parent;REPO=Path('/Users/yunhengz/Desktop/AAM Writing');M=V/'manuscript';OUT=V/'report'
OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
load=lambda p:json.loads(p.read_text())
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
d=load(V/'archive-metadata.json');manifest=d['manifest.json'];audit=V.parent/'current-validation/golden-data/audit.jsonl'
assert sha(audit)==manifest['audit_sha256']
assert sha(V/'archive-inputs.jsonl')==manifest['input_sha256']
current=[json.loads(x) for x in audit.read_text().splitlines()];inputs=[json.loads(x) for x in (V/'archive-inputs.jsonl').read_text().splitlines()]
assert inputs==[dict(index=r['index'],input_reaction=r['input_reaction']) for r in current]
assert len(current)==1851
original=load(REPO/'reports/golden_competitors_20260908/summary.json')
def function(path,name):return ast.dump(next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name))
assert function(V/'archived_evaluation_driver.py','signatures')==function(REPO/'bench/golden_competitors.py','signatures')
methods=[('rxnmapper','RXNMapper','0.4.3'),('localmapper','LocalMapper','0.1.5'),('chython','Chython','2.18 / chython-rxnmap 2.0'),('slap_binary','SLAP binary','1.0.0'),('slap_weighted','SLAP weighted','1.0.0'),('indigo','Indigo','1.46.0'),('rdt','Reaction Decoder Tool (RDT)','4.0.0')]
results=[];proofs={}
for method,name,version in methods:
 recheck=load(V/f'{method}-recheck.json');rows=recheck['records'];assert not recheck['changes'],recheck['changes'][:3]
 assert len(rows)==1851 and sorted(r['index'] for r in rows)==list(range(1851))
 first=[r['index'] for r in rows if r.get('first_correct')];any_=[r['index'] for r in rows if r.get('any_correct')]
 invalid=Counter(x['error'] for r in rows for x in r.get('invalid',[]));errors=[r['index'] for r in rows if r['status']!='mapped'];statuses=dict(Counter(r['status'] for r in rows))
 n=sum(r.get('candidates',0) for r in rows)
 assert len(first)==original[method]['first_correct'] and len(any_)==original[method]['any_correct']
 assert sum(invalid.values())==original[method]['invalid_candidates']
 assert statuses==original[method]['statuses']
 results.append(dict(method=method,name=name,version=version,total=1851,first_correct=len(first),any_correct=len(any_),first_correct_cases=first,any_correct_cases=any_,invalid_candidates=sum(invalid.values()),invalid_reasons=dict(invalid),failed_cases=errors,statuses=statuses,returned_candidates=n,original_mapping_timing={k:v for k,v in original[method].items() if 'seconds' in k or 'hours' in k}))
 proofs[method]={k:v for k,v in recheck.items() if k!='records'}
 with gzip.open(OUT/f'{method}-rescored.json.gz','wt') as f:json.dump(recheck,f,separators=(',',':'))
 package=dict(schema_version=1,denominator=1851,source='reports/golden_competitors_20260908/saved_outputs.tar.gz',source_sha256=sha(REPO/'reports/golden_competitors_20260908/saved_outputs.tar.gz'),original_summary_sha256=sha(REPO/'reports/golden_competitors_20260908/summary.json'),dataset_audit_sha256=sha(audit),unlabeled_inputs_sha256=manifest['input_sha256'],rescored_with_current_evaluator=True,mappers_rerun=False,unchanged_case_outcomes=True,methods=results,scope='Archived default-implementation outputs, freshly rescored. First and any explicitly returned mapping are separate. This is not exhaustive SLAP label-family decoding. All errors and invalid predictions stay in the denominator. No cross-host timing comparison is made.')
# Require every method before publishing the aggregate.
assert len(package['methods'])==7
save(M/'evidence/competitors.json',package);save(OUT/'competitors.json',package)
proof=dict(status='passed',checked_records=7*1851,workers=4,original_dataset_matches_current=True,signature_function_ast_unchanged=True,mapper_reruns=0,regression_tests='7 passed in 4.62s',evaluators={str(p.relative_to(REPO)):sha(p) for p in [REPO/'bench/golden_competitors.py',REPO/'bench/golden_evaluation.py']},driver_sha256=sha(V/'recheck_competitors.py'),method_rechecks=proofs)
save(OUT/'audit.json',proof)
for n in ['manifest.json','final_evaluation_provenance.json','competitor_env_20260908_provenance.json','neural_competitor_env_20260908_provenance.json']:save(OUT/n,d[n])
for n in ['recheck_competitors.py','assemble_comparators.py']:shutil.copy2(V/n,OUT/n)
(OUT/'README.md').write_text('''# Golden default-comparator recheck

All 12,957 saved method/reaction records from the seven completed default configurations were checked using the current strict heavy-atom relation evaluator. Mapping implementations and models were not rerun. Every per-case evaluation agrees with the archived result. The labeled audit and label-free input hashes match the final MAPPA/SLAP campaign exactly. Seven evaluator regression tests pass.

`competitors.json` supplies the paper table; compressed per-case files, `audit.json`, original version/model provenance, and the recheck drivers retain the evidence. Original predictions and failed attempts remain in `reports/golden_competitors_20260908/saved_outputs.tar.gz`, identified by SHA256. Signature evaluation preserves unmatched atoms and endpoint chemistry and handles agent-field spectators on the input side. Duplicate heavy labels are invalid.

The table distinguishes the first explicit mapping from recovery in any explicit returned mapping. It is not exhaustive SLAP symmetry-family decoding. MAPPA family recovery and the broader bidirectional SLAP sweep use separate search protocols. Original cluster mapping timings are retained in JSON but are not pooled with same-Mac MAPPA timing. RDT retries already completed in the original run; no final mapping failures remain for RDT.

Other archived tests are not conflated with these results. The coordinate-to-SMILES trials were cancelled and withdrawn; no scores from them are included. SAMMNet has no completed reproduction in this archive and receives no score. Superseded MAPPA experiments remain research history; the paper describes the final algorithm and current measured configurations.
''')
sources=load(M/'evidence/paper_sources.json');sources['snapshots']=[r for r in sources['snapshots'] if r['snapshot']!='competitors.json'];sources['snapshots'].append(dict(snapshot='competitors.json',source='reports/golden_competitor_recheck_20260913/competitors.json',scope=package['scope'],sha256=sha(M/'evidence/competitors.json')));save(M/'evidence/paper_sources.json',sources)
save(OUT/'SHA256SUMS.json',{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.json'})
print(json.dumps({'verified_records':7*1851,'methods':[(r['name'],r['first_correct'],r['any_correct']) for r in results]}))
