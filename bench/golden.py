"""Golden input preparation and reference policy, using the shared campaign runner."""
import argparse
from collections import defaultdict
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

from mappa import AAMSearchConfig, query_path
from mappa.mappa_core.core.scheduler.bundles import Blob, BundleWriter, read_blob
from mappa.mappa_core.utils.io.records import problem_from_record


REFERENCE_POLICY = 'exact heavy-atom relation including unmatched atoms, modulo endpoint chemical symmetry'


def endpoint_generators(feature, degree):
    """Lift chemically colored heavy-atom symmetries into original atom indexing."""
    import pynauty
    heavy = feature['heavy']
    colors, adjacency = defaultdict(set), defaultdict(set)
    for atom, color in enumerate(feature['colors']):
        colors[('atom', tuple(color))].add(atom)
    count = len(heavy)
    for left, right, color in feature['bonds']:
        colors[('bond', tuple(color))].add(count)
        adjacency[count].update((left, right))
        adjacency[left].add(count); adjacency[right].add(count)
        count += 1
    graph = pynauty.Graph(count, adjacency_dict={i: sorted(adjacency[i]) for i in range(count)},
        vertex_coloring=[colors[key] for key in sorted(colors, key=repr)])
    result = []
    for generator in pynauty.autgrp(graph)[0]:
        images = list(range(degree))
        for atom, original in enumerate(heavy):
            images[original] = heavy[generator[atom]]
        result.append(tuple(images))
    return tuple(result)


def verify(payload, workspace):
    saved = workspace.load(payload['catalogue'])
    reference = read_blob(payload['reference'], workspace.scratch)
    mapping = dict(reference['mapping'])
    features, generators = reference['features'], reference['generators']
    if payload['reversed']:
        mapping = {product: reactant for reactant, product in mapping.items()}
        features, generators = features[::-1], generators[::-1]
    unknown = False
    for number, family in enumerate(saved.catalogue.families):
        status, evidence = query_path(family.as_path(saved.catalogue.problem), saved.catalogue.problem,
            mapping, source_atoms=tuple(features[0]['heavy']), source_generators=generators[0],
            target_generators=generators[1], complete_reference=True, timeout_ms=10000)
        if status == 'recovered':
            return dict(status=status, family=number, evidence=evidence,
                        heavy_atoms=len(mapping), policy=REFERENCE_POLICY)
        unknown |= status == 'unknown'
    return dict(status='unknown' if unknown else 'not_recovered',
                heavy_atoms=len(mapping), policy=REFERENCE_POLICY)


def prepare(inputs, manifest, output, *, batches=5):
    """Copy verified original endpoints into aggregated bundles; never rerun old code."""
    inputs, manifest, output = Path(inputs), Path(manifest), Path(output)
    records = json.loads(manifest.read_text())['records']
    if sorted(row['case'] for row in records) != list(range(1851)):
        raise ValueError('Golden must contain exactly the original 1,851 records')
    output.mkdir(parents=True, exist_ok=False)
    writer = BundleWriter(output/'inputs', scratch=None)
    rows, evidence = [], []
    for row in sorted(records, key=lambda row: row['case']):
        folder = inputs/str(row['case'])
        raw = (folder/'input.json').read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != row['input_sha256']:
            raise ValueError(f'Golden input checksum differs for {row["case"]}')
        source = json.loads(raw)
        problem = problem_from_record(source)
        raw_reference = (folder/'reference.json').read_bytes()
        reference = json.loads(raw_reference)
        reference['generators'] = [endpoint_generators(feature, len(endpoint.elements))
            for feature, endpoint in zip(reference['features'], (problem.reactant, problem.product))]
        reference_blob = asdict(writer.put(reference))
        for reverse in (False, True):
            record = dict(source, reactant=source['product'], product=source['reactant']) if reverse else source
            rows.append(dict(case=row['case'], reversed=reverse,
                             input=asdict(writer.put(record)), reference=reference_blob))
        evidence.append(dict(case=row['case'], input_sha256=digest,
            reference_sha256=hashlib.sha256(raw_reference).hexdigest()))
    config = asdict(AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1.,
                                   random_seed=42, growth_policy='existence'))
    for batch in range(batches):
        folder = output/str(batch); folder.mkdir()
        subset = [row for row in rows if row['case'] % batches == batch]
        (folder/'index.json').write_text(json.dumps(dict(config=config, directions=subset,
            verification='golden:verify', reference_policy=REFERENCE_POLICY), indent=2)+'\n')
    (output/'manifest.json').write_text(json.dumps(dict(cases=1851, directions=len(rows), batches=batches,
        input_manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),
        config=config, reference_policy=REFERENCE_POLICY, records=evidence), indent=2)+'\n')
    print(json.dumps(dict(cases=1851, directions=len(rows), batches=batches, output=str(output))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs'); parser.add_argument('manifest'); parser.add_argument('output')
    parser.add_argument('--batches', type=int, default=5)
    args = parser.parse_args()
    prepare(args.inputs, args.manifest, args.output, batches=args.batches)
