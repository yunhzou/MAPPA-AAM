"""Prepare an isolated source copy; never modify production source."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--skip-build',action='store_true');args=p.parse_args()
package=Path(__file__).resolve().parent
work=Path(os.environ.get('MAPPA_EXPERIMENT_WORK',package/'work')).resolve()
repo=args.repo.resolve();manifest=json.loads((package/'manifest.json').read_text())
original=(repo/'src/mappa/competition.py').read_bytes()
if hashlib.sha256(original).hexdigest()!=manifest['original_competition_sha256']:
 raise SystemExit('Source changed: use the source_commit recorded in manifest.json.')
engine=work/'engine'
if engine.exists():raise SystemExit('Engine already exists; choose a new MAPPA_EXPERIMENT_WORK directory.')
engine.mkdir(parents=True)
shutil.copytree(repo/'src',engine/'src',ignore=shutil.ignore_patterns('__pycache__'))
for name in ['setup.py','pyproject.toml','README.md']:shutil.copy2(repo/name,engine/name)
(work/'competition-original.py').write_bytes(original)
subprocess.run(['patch','-p1','-i',str(package/'competition.patch')],cwd=engine,check=True)
assert hashlib.sha256((engine/'src/mappa/competition.py').read_bytes()).hexdigest()==manifest['experimental_competition_sha256']
if not args.skip_build:subprocess.run([sys.executable,'setup.py','build_ext','--inplace'],cwd=engine,check=True)
print(engine)
