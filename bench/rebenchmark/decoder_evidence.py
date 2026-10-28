"""Check saved event witnesses and completeness evidence without decoding again.

The caller must separately establish the expected task grid and saved-graph
catalogue identity. Solver limits here leave audit obligations unresolved.
"""
from collections import OrderedDict
from dataclasses import asdict
from functools import lru_cache
import hashlib
from pathlib import Path
import tempfile
import time

from mappa.event_patterns import SignedEventIndex
from mappa.final_branches import FinalFamily
from mappa.search_graph import frozen_value


class ProofReplay:
    """Reuse only conclusive checks of byte-identical SMT queries."""
    def __init__(self, seconds=10., save=lambda record: None, cache_size=4096):
        self.seconds, self.save, self.cache_size = seconds, save, cache_size
        self.cache = OrderedDict()

    def check(self, text, expected):
        import z3
        if expected not in ('sat', 'unsat'):
            raise ValueError('Only conclusive recorded queries can certify a result')
        digest = hashlib.sha256(text.encode()).hexdigest()
        if digest in self.cache:
            result = self.cache[digest]
            self.cache.move_to_end(digest)
            if result['status'] != expected:
                raise ValueError('Conflicting results for the same query')
            return result
        solver = z3.Solver()
        solver.from_string(text)
        solver.set(timeout=max(1, int(1000*self.seconds)))
        started = time.process_time()
        status = str(solver.check())
        result = dict(sha256=digest, status=status, expected=expected,
                      cpu_seconds=time.process_time()-started)
        if status == 'unknown':
            result['reason'] = solver.reason_unknown()
        elif status != expected:
            raise ValueError('Saved proof or witness does not replay')
        result['evidence'] = self.save(dict(result, smt2=text))
        if status != 'unknown':
            self.cache[digest] = result
            if len(self.cache) > self.cache_size:
                self.cache.popitem(last=False)
        return result


def validate_record(value, problem, read_query, replay):
    """Validate every stored candidate and each recorded family certificate.

    Candidate mappings are pinned exactly, without endpoint equivalence or
    reference mappings. Identical family representatives need no SMT query.
    """
    from mappa.family_query import compile_path
    decoded = value['decoded']
    config = decoded['config']
    index = SignedEventIndex(problem, threshold=config['threshold'],
                             metal_threshold=config['metal_threshold'])
    if frozen_value(decoded['event_policy']) != frozen_value(index.policy):
        raise ValueError('Event policy differs')
    families = [FinalFamily.from_record(row) for row in value['families']]
    reports = decoded['family_reports']
    if len(families) != len(reports) or len(families) != value['selected_families']:
        raise ValueError('Family/report coverage differs')
    source_edges = {}
    @lru_cache(maxsize=128)
    def path_for(number):
        family = families[number]
        floor = family.policy[0]
        if floor not in source_edges:
            source_edges[floor] = tuple((int(a),int(b)) for a,b in zip(index.a,index.b)
                                       if index.r[a,b] >= floor)
        return family.as_path(problem, source_edges[floor])
    proof_rows = []
    for number, (family, report) in enumerate(zip(families, reports)):
        family.validate_representative(problem)
        path = path_for(number)
        if report['family_id'] != number:
            raise ValueError('Family identity differs')
        window = decoded['event_window']
        bound = report['max_events']
        if bound is not None and (window is None or bound < window):
            raise ValueError('Family proof does not cover the saved window')
        references = report.get('query_evidence', [])
        if len(references) != report['solver_queries']:
            raise ValueError('Solver evidence is missing')
        statuses = []
        for reference in references:
            query = read_query(reference)
            if hashlib.sha256(query['smt2'].encode()).hexdigest() != query['sha256']:
                raise ValueError('Saved query checksum differs')
            statuses.append(query['status'])
            if query['status'] in ('sat', 'unsat'):
                proof_rows.append(replay.check(query['smt2'], query['status']))
            elif query['status'] != 'unknown':
                raise ValueError('Unexpected query status')
        if not report['complete']:
            continue
        reason = report['reason']
        vector = tuple(index.vector(dict(family.mapping)))
        certificates = index.family_certificates
        if reason == 'signed_event_invariance':
            valid = index.invariant_path(path)
        elif reason == 'source_target_class_invariance':
            valid = certificates.invariant_class(vector, certificates.actions(path))
        elif reason == 'reachable_image_event_lower_bound':
            valid = bound is not None and certificates.lower_bound(vector, certificates.actions(path), bound) > bound
        elif reason == 'all_event_orbits_excluded':
            valid = bool(statuses) and statuses[-1] == 'unsat' and all(s == 'sat' for s in statuses[:-1])
        else:
            valid = False
        if not valid:
            raise ValueError('Family completeness certificate failed')

    ids = set()
    for candidate in decoded['candidates']:
        mapping = {int(a): int(b) for a,b in candidate['mapping'].items()}
        if len(set(mapping.values())) != len(mapping):
            raise ValueError('Candidate reuses target atoms')
        if any(a < 0 or a >= problem.source_atom_count or b < 0 or b >= problem.target_atom_count
               or problem.reactant.elements[a] != problem.product.elements[b] for a,b in mapping.items()):
            raise ValueError('Invalid candidate atom correspondence')
        actual = index.describe(index.vector(mapping))
        for key in ('id', 'total', 'events', 'endpoint_events'):
            if frozen_value(actual[key]) != frozen_value(candidate[key]):
                raise ValueError('Candidate event witness differs: '+key)
        if candidate['id'] in ids:
            raise ValueError('Duplicate event candidate')
        ids.add(candidate['id'])
        if decoded['event_window'] is not None and candidate['total'] > decoded['event_window']:
            raise ValueError('Candidate lies outside the saved window')
        members = candidate['family_ids']
        if not members or any(i < 0 or i >= len(families) for i in members):
            raise ValueError('Candidate has no valid supporting family')
        family = families[members[0]]
        if set(mapping) != {a for a,_ in family.mapping}:
            raise ValueError('Candidate mapped source region differs')
        if mapping != dict(family.mapping):
            compiled = compile_path(path_for(members[0]), problem, mapping,
                                    source_atoms=tuple(sorted(mapping)), complete_reference=False)
            proof_rows.append(replay.check(compiled.solver.to_smt2(), 'sat'))

    minimum = min((c['total'] for c in decoded['candidates']), default=None)
    complete = decoded['incomplete_paths'] == 0 and all(r['complete'] for r in reports)
    proven = minimum is not None and decoded['incomplete_paths'] == 0 and (minimum == 0 or complete)
    if (decoded['observed_minimum'], decoded['complete'], decoded['minimum_proven']) != (minimum, complete, proven):
        raise ValueError('Saved minimum/completeness flags differ from their evidence')
    unknown = [r for r in proof_rows if r['status'] == 'unknown']
    return dict(verified=not unknown, candidates=len(ids), families=len(families),
                decoded_complete=complete, minimum=minimum, minimum_proven=proven,
                incomplete_families=[r['family_id'] for r in reports if not r['complete']],
                proof_checks=proof_rows, unresolved_proof_checks=len(unknown),
                scope='Saved-record evidence only; task coverage and catalogue/source identity must also be verified.')


def window_slack(value):
    """Both recorded decoder window policies have explicit, exclusive settings."""
    config = value['decoded']['config']
    if config['event_slack'] is not None:
        if config['max_events'] is not None or 'window_basis' in value:
            raise ValueError('Relative window contains conflicting absolute metadata')
        return config['event_slack']
    return value['window_basis']['slack']


def verify_window(value, problem, graph, tier):
    from .decoding import representative_bound
    config = value['decoded']['config']
    slack = window_slack(value)
    if config['event_slack'] is not None:
        minimum = value['decoded']['observed_minimum']
        bound = None if minimum is None else minimum + slack
        if value['decoded']['event_window'] != bound:
            raise ValueError('Relative minimum window differs')
    else:
        bound = representative_bound(problem, graph, value['seed']-1, tier, slack)
        if bound != value['window_basis']['bound'] or bound != config['max_events']:
            raise ValueError('Seed-prefix decoding bound differs')


def verify_task(value, task, expected_signature):
    """Use the immutable task, including metadata not repeated in its output."""
    from .decoder_inventory import signature
    key = f'{value["task_id"]:05d}/cut/{value["cut"]:05d}/seed/{value["seed"]-1:02d}'
    if (task.key != key or task.operation != 'rebenchmark.decoding:decode_seed' or
            signature(task) != expected_signature or value['search_rerun'] is not False):
        raise ValueError('Published output or expected task identity differs')
    actual = dict(source=value['source'], graph=value['graph'], raw=value['raw'],
        finalized=value['finalized'] is None, cut=value['cut'], task_id=value['task_id'],
        seed=value['seed']-1, slack=window_slack(value),
        family_seconds=value['decoded']['config']['seconds_per_family'])
    if any(frozen_value(task.payload[k]) != frozen_value(v) for k,v in actual.items()):
        raise ValueError('Output differs from the expected scientific task')
    return key


def audit_output(payload, workspace):
    """Shared-executor operation for one accepted saved output, with provenance."""
    from mappa.artifacts import read_graph_checkpoint
    from mappa.execution import Blob, Task, read_blob
    from mappa.final_branches import FinalBranchCatalogue
    from .campaign import load_plan
    from .decoding import coverage
    from .saved_analysis import read_member
    from .search import select_contexts

    blob = Blob(**payload['blob'])
    value = read_blob(blob, workspace.scratch)
    key = verify_task(value, Task(**payload['specification']), payload['specification_signature'])
    if key != payload['key']:
        raise ValueError('Published output identity differs')
    config = value['decoded']['config']

    def prepare():
        with tempfile.TemporaryDirectory(dir=workspace.scratch) as td:
            folder = Path(td)
            for name, member in value['source'].items():
                (folder/name).write_bytes(read_member(member))
            plan = load_plan(folder)
            read_member(value['raw'])  # Verify the immutable raw-cut source too.
            if value['finalized'] is None:
                data = read_member(value['graph'])
                if hashlib.sha256(data).hexdigest() != value['finalized_sha256']:
                    raise ValueError('Finalized graph checksum differs')
                path = folder/'graph.pkl.gz'; path.write_bytes(data)
                graph = read_graph_checkpoint(path)
            else:
                final = value['finalized']
                if final.sha256 != value['finalized_sha256']:
                    raise ValueError('Saved finalized graph identity differs')
                schema, graph = read_blob(final, workspace.scratch)
                if schema != 'rxn_core.finalized_cut/v1':
                    raise ValueError('Unexpected finalized graph schema')
            return plan, graph
    plan, graph = workspace.prepared(
        ('audit-cut', value['source']['task.json']['sha256'], value['raw']['sha256'],
         value['finalized_sha256']), prepare)
    catalogue = FinalBranchCatalogue(plan.problem).add_graph(select_contexts(graph,[value['seed']-1]))
    tier = max((coverage(plan.problem,f) for f in catalogue.families), default=None)
    families = [f.to_record() for f in catalogue.families if coverage(plan.problem,f) == tier]
    if (len(catalogue.families) != value['all_families'] or
            frozen_value(families) != frozen_value(value['families']) or
            frozen_value(tier) != frozen_value(value['coverage']) or
            catalogue.incomplete_paths != value['decoded']['incomplete_paths']):
        raise ValueError('Saved family catalogue differs from the saved search graph')
    verify_window(value, plan.problem, graph, tier)
    if not hasattr(workspace, 'proof_replay'):
        workspace.proof_replay = ProofReplay(seconds=payload['proof_seconds'],
            save=lambda record: asdict(workspace.writer.put(record)))
    result = validate_record(value, plan.problem,
        lambda ref: read_blob(Blob(**ref), workspace.scratch), workspace.proof_replay)
    result.update(key=key, blob=payload['blob'], specification_signature=payload['specification_signature'],
                  catalogue_verified=True, prefix_bound_verified=True,
                  scope='This accepted output only; dataset coverage requires every expected key.')
    return result
