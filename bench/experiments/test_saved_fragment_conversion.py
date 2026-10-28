"""One saved seed/cut conversion; no search, solver, or shared extraction."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import tarfile
import tempfile

from mappa.artifacts import read_graph_checkpoint
from mappa.final_branches import FinalBranchCatalogue
from mappa.fragment_conversion import convert_aam_fragments
from rebenchmark.campaign import load_plan
from rebenchmark.search import select_contexts


def run(archive, cut, seed, output):
    member = f'finalized/{cut:05d}.pkl.gz'
    with tempfile.TemporaryDirectory() as scratch, tarfile.open(archive) as saved:
        folder = Path(scratch)
        manifest = json.load(saved.extractfile('archive-manifest.json'))['files']
        hashes = {}
        for name in ('input.json', 'task.json', member):
            data = saved.extractfile(name).read()
            digest = hashlib.sha256(data).hexdigest()
            assert digest == manifest[name]['sha256'] and len(data) == manifest[name]['bytes']
            path = folder / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            hashes[name] = digest
        plan = load_plan(folder)
        graph = select_contexts(read_graph_checkpoint(folder / member), [seed])
        catalogue = FinalBranchCatalogue(plan.problem).add_graph(graph, str(archive))
        before = catalogue.to_record()
        result = convert_aam_fragments(catalogue)
        assert before == catalogue.to_record() == result.to_record()['original_catalogue']
        # Independent dense bond-difference check, including partial boundaries.
        verified = 0
        for branch in result.branches:
            if not branch.complete:
                continue
            mapping = dict(catalogue.families[branch.family_id].mapping)
            problem = catalogue.problem
            from mappa.frag import bond_event_threshold
            events = 0
            for a in range(problem.source_atom_count):
                for b in range(a + 1, problem.source_atom_count):
                    if a not in mapping and b not in mapping:
                        continue
                    r = problem.reactant.wbo[a, b]
                    p = problem.product.wbo[mapping[a], mapping[b]] if a in mapping and b in mapping else 0
                    tau = bond_event_threshold(problem.reactant.elements, a, b,
                                               default_threshold=.5, metal_threshold=.3)
                    events += abs(p - r) >= tau
            targets = set(mapping.values())
            for a in range(problem.target_atom_count):
                for b in range(a + 1, problem.target_atom_count):
                    if (a in targets) == (b in targets):
                        continue
                    tau = bond_event_threshold(problem.product.elements, a, b,
                                               default_threshold=.5, metal_threshold=.3)
                    events += problem.product.wbo[a, b] >= tau
            assert branch.events == events
            verified += 1
        record = result.to_record()
        record.update(input_archive=str(archive), input_hashes=hashes, cut=cut, seed_index=seed,
                      job_id=os.environ.get('SLURM_JOB_ID'),
                      independent_counts_verified=verified,
                      unresolved_reasons=dict(Counter(b.reason for b in result.branches if not b.complete)))
        output.write_text(json.dumps(record, indent=2) + '\n')
        print(json.dumps({k: v for k, v in record.items() if k not in ('original_catalogue', 'branches')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--cut', type=int, default=0)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.archive, args.cut, args.seed, args.output)
