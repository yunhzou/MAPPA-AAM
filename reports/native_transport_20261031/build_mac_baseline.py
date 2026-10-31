"""Build the pre-refactor transport for the local point test without changing the checkout."""
import importlib.util
from pathlib import Path
import subprocess
import sys

repo=Path(__file__).resolve().parents[2]
before=Path(sys.argv[1]).resolve()
before.mkdir(parents=True,exist_ok=True)
files={
    'engine.cpp':'src/mappa/mappa_core/native/engine.cpp',
    'fragment_scheduler.h':'src/mappa/mappa_core/native/fragment_scheduler.h',
    'graph.py':'src/mappa/mappa_core/core/search/graph.py',
    'adapter.py':'src/mappa/mappa_core/core/search/branching/native.py',
}
for name,path in files.items():
    (before/name).write_bytes(subprocess.check_output(['git','show','733adfc:'+path],cwd=repo))
spec=importlib.util.spec_from_file_location('build_engine',repo/'src/mappa/mappa_core/native/build_engine.py')
build=importlib.util.module_from_spec(spec); spec.loader.exec_module(build)
build.ext.sources[0]=str(before/'engine.cpp')
build.ext.include_dirs.append(str(repo/'src/mappa/mappa_core/native'))
build.ROOT=before
(before/'src/mappa').mkdir(parents=True,exist_ok=True)
build.main()
