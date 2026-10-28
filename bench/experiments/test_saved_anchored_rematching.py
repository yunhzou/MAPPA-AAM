"""Rematch saved fragment pairs and compare recombined branch event counts."""
import argparse
from collections import Counter
import gzip
import hashlib
import io
from itertools import product
import json
import os
from pathlib import Path
import tarfile
import time

from mappa import AAMProblem, MolecularEndpoint, BondTolerance
from mappa.alignment.branch import find_islands
from mappa.artifacts import _load_checkpoint
from mappa.event_patterns import SignedEventIndex
from mappa.final_branches import FinalBranchCatalogue
from mappa.frag import bond_event_threshold, build_graph, classify_bonds
from mappa.search_symmetry import finalize_graph_symmetry


def run(archive, expected_sha, output):
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == expected_sha
    with tarfile.open(archive) as saved:
        manifest = json.load(saved.extractfile('__mappa_archive_manifest__.json'))['files']
        def read(name):
            data = saved.extractfile(name).read()
            assert hashlib.sha256(data).hexdigest() == manifest[name]['sha256']
            return data
        raw = json.loads(gzip.decompress(read('input.json.gz')))
        problem = AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant', 'product')), raw['name'])
        schema, fields = _load_checkpoint(io.BytesIO(gzip.decompress(read('catalogue.shared.pkl.gz'))))
        assert schema == 'mappa.shared_catalogue/v1'
        catalogue = FinalBranchCatalogue(problem)
        assert fields['problem_sha256'] == catalogue.problem_sha256
        vars(catalogue).update(fields)
        prior = json.loads(gzip.decompress(read('decoding/definitive-decode-v1/result.json.gz')))
        assert prior['complete'] and prior['minimum_proven']

    threshold, metal = prior['config']['threshold'], prior['config']['metal_threshold']
    def tolerance(graph, a, b):
        return bond_event_threshold(problem.reactant.elements, a, b,
                                    default_threshold=threshold, metal_threshold=metal)
    rule = BondTolerance(tolerance, inclusive=False)
    graphs, cache, pair_rows = {}, {}, []
    entries = [(policy, pairs, fid) for (policy, pairs), ids in catalogue.branches.items()
               for fid in sorted(ids)]
    entries += [(catalogue.families[fid].policy, pairs, fid)
                for fid, pairs_list in catalogue.coupled_branches.items() for pairs in sorted(pairs_list)]
    started, cpu = time.perf_counter(), time.process_time()
    histogram, best, rows = Counter(), [], []
    calls = reused = internal_events = incomplete = 0
    for policy, pairs, fid in entries:
        family = catalogue.families[fid]
        anchors = dict(family.fixed_assignments)
        if policy[0] not in graphs:
            graphs[policy[0]] = tuple(build_graph(ep.elements, ep.wbo, bond_cut=policy[0])
                                    for ep in (problem.reactant, problem.product))
        choices = []
        for rs, ps in pairs:
            fixed = tuple((r, anchors[r]) for r in rs if r in anchors)
            key = (policy[0], rs, ps, fixed)
            if key not in cache:
                source, target = (g.subgraph(atoms).copy() for g, atoms in zip(graphs[policy[0]], (rs, ps)))
                pair_start = time.perf_counter()
                result = find_islands(source, target, sorted(rs),
                    graph_floor=policy[0], iso_tol=threshold, max_branches=100,
                    anchor_map=dict(fixed), bond_tolerance=rule, random_seed=42)
                result, _ = finalize_graph_symmetry(result, target, iso_tolerance=threshold)
                paths = list(result.paths())
                states = tuple(result.states[p.terminal] for p in paths)
                complete = not result.capped and bool(states) and all(len(s.mapping) == len(rs) for s in states)
                cache[key] = states
                pair_rows.append(dict(source_atoms=rs, target_atoms=ps, anchors=fixed,
                    branches=len(states), capped=result.capped, complete=complete,
                    seconds=time.perf_counter()-pair_start, search=result.to_record()))
                incomplete += not complete
                calls += 1
            else:
                reused += 1
            choices.append(cache[key])
        for states in product(*choices):
            mapping = {r: p for s in states for r, p in s.mapping}
            assert len(mapping) == sum(len(s.mapping) for s in states)
            assert len(set(mapping.values())) == len(mapping)
            assert all(mapping[r] == p for r, p in anchors.items())
            broken, formed, _, _ = classify_bonds(mapping, problem.reactant.wbo, problem.product.wbo,
                dwbo_threshold=threshold, metal_dwbo_threshold=metal,
                elements_R=problem.reactant.elements, elements_P=problem.product.elements)
            # Independent dense subtraction, including mapped/unmapped boundaries.
            direct = []
            for left, right, assignment in ((problem.reactant, problem.product, mapping),
                                            (problem.product, problem.reactant, {p:r for r,p in mapping.items()})):
                events = set()
                for a in range(left.atom_count):
                    for b in range(a+1, left.atom_count):
                        if a not in assignment and b not in assignment:
                            continue
                        other = right.wbo[assignment[a], assignment[b]] if a in assignment and b in assignment else 0.
                        limit = bond_event_threshold(left.elements, a, b, default_threshold=threshold, metal_threshold=metal)
                        if left.wbo[a,b] - other >= limit:
                            events.add((a,b))
                direct.append(events)
            assert direct == [set((a,b) for a,b,*_ in broken), set((a,b) for a,b,*_ in formed)]
            owners = {r:(i,island) for i,s in enumerate(states) for r,island in s.islands}
            target_owners = {mapping[r]:owner for r,owner in owners.items()}
            internal = sum(a in own and b in own and own[a] == own[b]
                           for events, own in ((broken, owners), (formed, target_owners)) for a,b,*_ in events)
            internal_events += internal
            count = len(broken)+len(formed)
            histogram[count] += 1
            row = dict(family=fid, events=count, mapped_atoms=len(mapping),
                       changed_assignments=sum(mapping.get(r) != p for r,p in family.mapping),
                       internal_events=internal)
            rows.append(row)
            if not best or count < best[0]['events']:
                best = []
            if not best or count == best[0]['events']:
                best.append(dict(row, mapping=mapping, broken=sorted(direct[0]), formed=sorted(direct[1])))
    rematch_seconds, rematch_cpu = time.perf_counter()-started, time.process_time()-cpu
    index = SignedEventIndex(problem, threshold=threshold, metal_threshold=metal)
    old_classes = {index.describe(index.vector({int(a):int(b) for a,b in c['mapping'].items()}))['id']
                   for c in prior['candidates'] if c['total'] == prior['observed_minimum']}
    descriptions = {}
    for row in best:
        vector = tuple(index.vector(row['mapping']))
        if vector not in descriptions:
            descriptions[vector] = index.describe(vector)['id']
        row['event_class'] = descriptions[vector]
    summary = dict(case=problem.name, source_atoms=problem.source_atom_count,
        threshold=threshold, metal_threshold=metal, inclusive=False, branch_cap=100,
        saved_families=len(catalogue.families), original_branch_entries=len(entries),
        distinct_fragment_pairs=len({(k[1],k[2]) for k in cache}), distinct_anchored_rematches=calls,
        reused_pair_results=reused, new_branch_combinations=len(rows), event_histogram=dict(sorted(histogram.items())),
        minimum_new_branch_events=min(histogram, default=None), prior_minimum=prior['observed_minimum'],
        best_branches_matching_prior_minimum_event_class=sum(r['event_class'] in old_classes for r in best),
        branches_with_changed_free_assignments=sum(r['changed_assignments'] > 0 for r in rows),
        all_branches_map_all_source_atoms=all(r['mapped_atoms'] == problem.source_atom_count for r in rows),
        incomplete_pair_rematches=incomplete, internal_events=internal_events,
        wall_seconds_including_independent_checks=rematch_seconds,
        cpu_seconds_including_independent_checks=rematch_cpu,
        pair_matching_seconds=sum(r['seconds'] for r in pair_rows), solver_queries=0,
        job_id=os.environ.get('SLURM_JOB_ID'), archive_sha256=expected_sha,
        scope='Direct counts across recombined new compressed-branch witnesses; not exhaustive old-decoder window equivalence.')
    output.write_text(json.dumps(dict(summary=summary, examples=best[:2], branches=rows), indent=2)+'\n')
    with gzip.open(output.with_suffix('.pairs.json.gz'), 'wt') as stream:
        json.dump(pair_rows, stream)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--archive-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.archive, args.archive_sha256, args.output)
