"""Compare representative refinement with an audited saved decoding result."""
import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile

from mappa import AAMProblem, MolecularEndpoint
from mappa.artifacts import _load_checkpoint
from mappa.event_patterns import SignedEventIndex
from mappa.final_branches import FinalBranchCatalogue
from mappa.fragment_conversion import convert_aam_fragments
from mappa.frag import bond_event_threshold


def direct_events(problem, mapping, threshold, metal_threshold):
    """Independent dense subtraction; return bonds in endpoint-native indices."""
    inverse = {p: r for r, p in mapping.items()}
    result = []
    for source, target, assignment in ((problem.reactant, problem.product, mapping),
                                       (problem.product, problem.reactant, inverse)):
        events = []
        for a in range(source.atom_count):
            for b in range(a + 1, source.atom_count):
                if a not in assignment and b not in assignment:
                    continue
                other = target.wbo[assignment[a], assignment[b]] if a in assignment and b in assignment else 0.
                limit = bond_event_threshold(source.elements, a, b,
                    default_threshold=threshold, metal_threshold=metal_threshold)
                if source.wbo[a, b] - other >= limit:
                    events.append((a, b))
        result.append(events)
    return tuple(result)


def compare(archive, archive_sha256, output):
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    assert digest == archive_sha256, 'Archive differs from the independent audit receipt'
    with tarfile.open(archive) as saved:
        manifest = json.load(saved.extractfile('__mappa_archive_manifest__.json'))['files']
        hashes = {}

        def member(name):
            data = saved.extractfile(name).read()
            hashes[name] = hashlib.sha256(data).hexdigest()
            assert hashes[name] == manifest[name]['sha256'] and len(data) == manifest[name]['size']
            return data

        raw = json.loads(gzip.decompress(member('input.json.gz')))
        problem = AAMProblem(*(MolecularEndpoint(**raw[side]) for side in ('reactant', 'product')), raw['name'])
        schema, fields = _load_checkpoint(io.BytesIO(gzip.decompress(member('catalogue.shared.pkl.gz'))))
        assert schema == 'mappa.shared_catalogue/v1'
        catalogue = FinalBranchCatalogue(problem)
        assert fields['problem_sha256'] == catalogue.problem_sha256
        assert set(fields) == set(vars(catalogue)) - {'problem'}
        vars(catalogue).update(fields)
        before = catalogue.to_record()
        prior = json.loads(gzip.decompress(member('decoding/definitive-decode-v1/result.json.gz')))
        source = json.loads(member('catalogue-progress.json'))
        assert source['complete'] and prior['complete'] and prior['minimum_proven']
        threshold, metal = prior['config']['threshold'], prior['config']['metal_threshold']
        converted = convert_aam_fragments(catalogue, threshold=threshold, metal_threshold=metal)
        assert catalogue.to_record() == before
        index = SignedEventIndex(problem, threshold=threshold, metal_threshold=metal)
        rows = []
        for branch in converted.branches:
            mapping = dict(catalogue.families[branch.family_id].mapping)
            broken, formed = direct_events(problem, mapping, threshold, metal)
            assert broken == [(a, b) for a, b, *_ in branch.broken]
            assert formed == [(a, b) for a, b, *_ in branch.formed]
            description = index.describe(index.vector(mapping))
            assert description['total'] == branch.representative_events
            rows.append(dict(family_id=branch.family_id, mapping=mapping,
                original_fragments=len(branch.original_pairs), refined_fragments=len(branch.fragment_pairs),
                fragment_pairs=branch.fragment_pairs, partition_complete=branch.partition_complete,
                broken=broken, formed=formed, events=branch.representative_events,
                event_class=description['id']))
        lowest = min(row['events'] for row in rows)
        best = [row for row in rows if row['events'] == lowest]
        old_classes = {c['id'] for c in prior['candidates'] if c['total'] == prior['observed_minimum']}
        summary = dict(case=problem.name, source_atoms=problem.source_atom_count,
            target_atoms=problem.target_atom_count, cuts=source['total_cuts'],
            saved_fragment_branches=catalogue.branch_count, saved_mapping_families=len(catalogue.families),
            threshold=threshold, metal_threshold=metal, prior_minimum=prior['observed_minimum'],
            prior_candidates=len(prior['candidates']), lowest_converted_representative_count=lowest,
            equal_minimum_count=lowest == prior['observed_minimum'],
            converted_count_histogram=dict(sorted(Counter(r['events'] for r in rows).items())),
            representatives_at_lowest_count=len(best),
            best_event_classes_matching_old=sum(r['event_class'] in old_classes for r in best),
            old_minimum_event_classes=len(old_classes),
            independently_checked_representatives=len(rows), statistics=converted.statistics,
            unresolved_partition_count=sum(not r['partition_complete'] for r in rows),
            scope='Original saved representatives refined at event thresholds; no expansion of free compressed choices.',
            complete_compressed_conversion=converted.complete,
            archive=str(archive), archive_sha256=digest, member_hashes=hashes,
            job_id=os.environ.get('SLURM_JOB_ID'))
        output.write_text(json.dumps(dict(summary=summary, example=best[0], branches=rows), indent=2)+'\n')
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--archive-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    compare(args.archive, args.archive_sha256, args.output)
