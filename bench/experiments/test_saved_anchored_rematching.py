"""Rematch saved fragment pairs and compare recombined branch event counts."""
import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import time

from mappa import AAMProblem, MolecularEndpoint
from mappa.artifacts import _load_checkpoint
from mappa.event_patterns import SignedEventIndex
from mappa.final_branches import FinalBranchCatalogue
from mappa.frag import bond_event_threshold, classify_bonds


def run(archive, expected_sha, output,
        result_member='decoding/definitive-decode-v1/result.json.gz'):
    digest = hashlib.sha256()
    with archive.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    assert digest.hexdigest() == expected_sha
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
        prior = json.loads(gzip.decompress(read(result_member)))
        assert prior['complete'] and prior['minimum_proven']

    from mappa.postprocessing import decode_events, EventDecodeConfig
    threshold, metal = prior['config']['threshold'], prior['config']['metal_threshold']
    started, cpu = time.perf_counter(), time.process_time()
    decoded = decode_events(catalogue, EventDecodeConfig(threshold=threshold, metal_threshold=metal))
    conversion = decoded.conversion
    conversion_seconds = time.perf_counter()-started
    pair_rows = []
    serialize_start = time.perf_counter()
    for pair in conversion.rematches:
        pair_rows.append(dict(source_atoms=pair.source_atoms, target_atoms=pair.target_atoms,
            anchors=pair.anchors, branches=len(pair.paths), capped=pair.graph.capped,
            complete=pair.complete, seconds=pair.seconds, search=pair.graph.to_record()))
    timings = dict(graph_serialization_seconds=time.perf_counter()-serialize_start,
                   conversion_seconds=conversion_seconds,
                   matching_seconds=sum(p.seconds for p in conversion.rematches),
                   rematching_recombination_seconds=conversion.seconds,
                   classification_and_grouping_seconds=conversion_seconds-conversion.seconds)
    # This intentionally separate dense oracle checks every distinct result map.
    check_start = time.perf_counter()
    check_pairs = tuple(tuple((a, b, ep.wbo[a, b],
        bond_event_threshold(ep.elements, a, b, default_threshold=threshold, metal_threshold=metal))
        for a in range(ep.atom_count) for b in range(a+1, ep.atom_count))
        for ep in (problem.reactant, problem.product))
    verified, histogram, rows, best = {}, Counter(), [], []
    reports = {r['family_id']:r for r in decoded.family_reports}
    candidates = {fid:c for c in decoded.candidates for fid in c.family_ids}
    for branch in conversion.branches:
        family = conversion.catalogue.families[branch.family_id]
        mapping = dict(family.mapping)
        candidate = candidates[branch.family_id]
        if family.mapping not in verified:
            direct = []
            for pairs_to_check, right, assignment in ((check_pairs[0], problem.product, mapping),
                                           (check_pairs[1], problem.reactant, {p:r for r,p in mapping.items()})):
                events = set()
                for a, b, weight, limit in pairs_to_check:
                    if a not in assignment and b not in assignment:
                        continue
                    other = right.wbo[assignment[a], assignment[b]] if a in assignment and b in assignment else 0.
                    if weight - other >= limit:
                        events.add((a,b))
                direct.append(events)
            # Candidate event identities may group different literal maps;
            # compare this branch's literal count and its own physical bonds.
            broken, formed, _, _ = classify_bonds(mapping, problem.reactant.wbo, problem.product.wbo,
                dwbo_threshold=threshold, metal_dwbo_threshold=metal,
                elements_R=problem.reactant.elements, elements_P=problem.product.elements)
            assert direct == [set((a,b) for a,b,*_ in broken), set((a,b) for a,b,*_ in formed)]
            assert sum(map(len,direct)) == candidate.total
            verified[family.mapping] = direct
        direct = verified[family.mapping]
        count = candidate.total
        histogram[count] += 1
        old = catalogue.families[branch.source_family_id]
        row = dict(family=branch.source_family_id, events=count, mapped_atoms=len(mapping),
                   changed_assignments=sum(mapping.get(r) != p for r,p in old.mapping),
                   internal_events=reports[branch.family_id]['internal_events'])
        rows.append(row)
        if not best or count < best[0]['events']:
            best = []
        if not best or count == best[0]['events']:
            best.append(dict(row, mapping=mapping, broken=sorted(direct[0]), formed=sorted(direct[1])))
    timings['independent_check_seconds'] = time.perf_counter()-check_start
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
        saved_families=len(catalogue.families), original_branch_entries=sum(len(ids) for ids in catalogue.branches.values()) + sum(len(pairs) for pairs in catalogue.coupled_branches.values()),
        distinct_fragment_pairs=len({(p.source_atoms,p.target_atoms) for p in conversion.rematches}), distinct_anchored_rematches=len(conversion.rematches),
        reused_pair_results=conversion.reused_pairs, unique_counted_mappings=len(verified),
        reused_mapping_counts=len(rows)-len(verified), new_branch_combinations=len(rows), event_histogram=dict(sorted(histogram.items())),
        minimum_new_branch_events=min(histogram, default=None), prior_minimum=prior['observed_minimum'],
        best_branches_matching_prior_minimum_event_class=sum(r['event_class'] in old_classes for r in best),
        branches_with_changed_free_assignments=sum(r['changed_assignments'] > 0 for r in rows),
        all_branches_map_all_source_atoms=all(r['mapped_atoms'] == problem.source_atom_count for r in rows),
        incomplete_pair_rematches=sum(not p.complete for p in conversion.rematches), internal_events=sum(r["internal_events"] for r in rows),
        wall_seconds_including_independent_checks=rematch_seconds,
        timings=dict(timings),
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
    parser.add_argument('--result-member', default='decoding/definitive-decode-v1/result.json.gz')
    args = parser.parse_args()
    run(args.archive, args.archive_sha256, args.output, args.result_member)
