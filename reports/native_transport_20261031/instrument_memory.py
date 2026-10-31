"""Build current native code with cache measurements and a diagnostic retention control.

The output stays in the requested diagnostic directory; no installed engine is replaced.
"""
import argparse
from pathlib import Path
import runpy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    args = parser.parse_args()
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=True)
    native = Path(__file__).resolve().parents[2]/'src/mappa/mappa_core/native'
    source = (native/'engine.cpp').read_text()
    marker = '        .def("memo_stats", [](const PyTarget& t) {'
    assert source.count(marker) == 1
    statistics = '''        .def("diagnostic_memo_limit", [](PyTarget& t, size_t entries) {
            t.g.nauty->memo_limit=entries;
        })
        .def("memory_stats", [](const PyTarget& t) {
            const auto& graph=*t.g.nauty;
            size_t keys=0, values=0, unique_values=0;
            std::unordered_set<std::string_view> distinct;
            for (const auto& item:graph.memo) {
                keys+=item.first.capacity()+1;
                values+=item.second.capacity()+1;
                if (distinct.emplace(item.second).second) unique_values+=item.second.size();
            }
            py::dict out;
            out["vertices"]=graph.n;
            out["certificate_bytes"]=graph.canon.size()*sizeof(setword);
            out["memo_entries"]=graph.memo.size();
            out["memo_key_bytes"]=keys;
            out["memo_value_bytes"]=values;
            out["unique_certificate_values"]=distinct.size();
            out["unique_certificate_bytes"]=unique_values;
            out["roles"]=t.g.roles.ids.size();
            out["colours"]=t.g.colours.ids.size();
            return out;
        })
'''
    source = '#include <string_view>\n'+source.replace(marker, statistics+marker)
    instrumented = output/'engine.cpp'; instrumented.write_text(source)
    build = runpy.run_path(str(native/'build_engine.py'))
    extension = build['ext']
    extension.sources[0] = str(instrumented)
    extension.include_dirs.append(str(native))
    distribution = build['Distribution']({'name':'mappa_memory_probe', 'ext_modules':[extension]})
    command = build['Build'](distribution)
    command.build_lib = str(output/'lib'); command.build_temp = str(output/'obj')
    command.ensure_finalized(); command.run()
    print(next((output/'lib/mappa').glob('_engine*.so')))


if __name__ == '__main__':
    main()
