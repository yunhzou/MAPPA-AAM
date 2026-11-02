"""Build current native code with one memory-ownership report before Python export.

The diagnostic engine stays in the requested directory. Search, returned records,
and the installed engine remain unchanged. Counts describe owned vector capacity;
allocator metadata, shared-pointer control blocks and string allocations are excluded.
"""
import argparse
from pathlib import Path
import runpy


HISTORY_MEASUREMENTS = r'''
    auto bytes=[](const auto& values) { return values.capacity()*sizeof(values[0]); };
    size_t state_mapping_labels=0,state_deferred=0,state_deferred_pairs=0;
    for (const auto& state:graph.states) {
        state_mapping_labels+=bytes(state.mapping)+bytes(state.labels);
        state_deferred+=bytes(state.deferred);
        state_deferred_pairs+=state.deferred.size();
    }
    size_t matches=0,candidate_images=0,match_deferred=0,match_deferred_pairs=0;
    size_t match_other_vectors=0,pools=0,transition_preserved=0;
    for (const auto& transition:graph.transitions) {
        transition_preserved+=bytes(transition.preserved);
        if (!transition.match) continue;
        // Each committed fragment owns a distinct IsoOut in this history.
        ++matches;
        const auto& match=*transition.match;
        candidate_images+=bytes(match.cand.img);
        match_deferred+=bytes(match.deferred_edges);
        match_deferred_pairs+=match.deferred_edges.size();
        match_other_vectors+=bytes(match.mapping)+bytes(match.fragment)+bytes(match.preserved_bonds)
            +bytes(match.cand.mapped)+bytes(match.cand.exact_fixed);
        pools+=bytes(match.cand.blocks);
        for (const auto& pool:match.cand.blocks) pools+=bytes(pool.r)+bytes(pool.p);
    }
    const size_t graph_storage=sizeof(graph)+bytes(graph.states)+bytes(graph.transitions)+bytes(graph.stops);
    const size_t match_objects=matches*sizeof(IsoOut);
    const size_t total=graph_storage+match_objects+state_mapping_labels+state_deferred
        +candidate_images+match_deferred+match_other_vectors+pools+transition_preserved;
    std::fprintf(stderr,
        "{\"phase\":\"native_history\",\"source_atoms\":%zu,\"target_atoms\":%zu,"
        "\"states\":%zu,\"transitions\":%zu,\"matches\":%zu,\"stops\":%zu,"
        "\"state_deferred_pairs\":%zu,\"match_deferred_pairs\":%zu,"
        "\"state_mapping_labels_bytes\":%zu,\"state_deferred_bytes\":%zu,"
        "\"candidate_image_bytes\":%zu,\"match_deferred_bytes\":%zu,"
        "\"match_other_vector_bytes\":%zu,\"pool_storage_bytes\":%zu,"
        "\"transition_preserved_bytes\":%zu,\"graph_object_and_vector_bytes\":%zu,"
        "\"match_object_bytes\":%zu,\"accounted_capacity_bytes\":%zu}\n",
        r_ids.size(),p_ids.size(),graph.states.size(),graph.transitions.size(),matches,graph.stops.size(),
        state_deferred_pairs,match_deferred_pairs,state_mapping_labels,state_deferred,
        candidate_images,match_deferred,match_other_vectors,pools,transition_preserved,
        graph_storage,match_objects,total);
    std::fflush(stderr);
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    native = Path(__file__).resolve().parents[2]/'src/mappa/mappa_core/native'
    binding = (native/'fragment_scheduler_binding.h').read_text()
    marker = '    SharedStatePairs shared;\n'
    assert binding.count(marker) == 1
    (output/'fragment_scheduler_binding.h').write_text(binding.replace(marker, HISTORY_MEASUREMENTS+marker))
    instrumented = output/'engine.cpp'
    instrumented.write_text('#include <cstdio>\n'+(native/'engine.cpp').read_text())
    build = runpy.run_path(str(native/'build_engine.py'))
    extension = build['ext']
    extension.sources[0] = str(instrumented)
    extension.include_dirs.append(str(native))
    distribution = build['Distribution']({'name': 'mappa_memory_probe', 'ext_modules': [extension]})
    command = build['Build'](distribution)
    command.build_lib = str(output/'lib')
    command.build_temp = str(output/'obj')
    command.ensure_finalized()
    command.run()
    print(next((output/'lib/mappa').glob('_engine*.so')))


if __name__ == '__main__':
    main()
