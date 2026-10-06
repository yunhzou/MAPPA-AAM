"""Build the prefix experiment in an isolated engine directory."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--repo',required=True,type=Path);p.add_argument('--skip-build',action='store_true');args=p.parse_args()
package=Path(__file__).resolve().parent;repo=args.repo.resolve()
work=Path(os.environ.get('MAPPA_EXPERIMENT_WORK',package/'work')).resolve();engine=work/'engine'
manifest=json.loads((package/'manifest.json').read_text())
for path,hashes in manifest['source_hashes'].items():
 if hashlib.sha256((repo/path).read_bytes()).hexdigest()!=hashes['original']:
  raise SystemExit('Source changed: use source_commit from manifest.json.')
if engine.exists():raise SystemExit('Choose a fresh MAPPA_EXPERIMENT_WORK directory.')
engine.mkdir(parents=True)
for name in ['src','native']:shutil.copytree(repo/name,engine/name,ignore=shutil.ignore_patterns('__pycache__'))
for name in ['setup.py','pyproject.toml','README.md']:shutil.copy2(repo/name,engine/name)
shutil.copy2(repo/'src/mappa/competition.py',work/'competition-original.py')
subprocess.run(['patch','-p1','-i',str(package/'prefix.patch')],cwd=engine,check=True)
for path,hashes in manifest['source_hashes'].items():
 assert hashlib.sha256((engine/path).read_bytes()).hexdigest()==hashes['experimental']
if not args.skip_build:
 subprocess.run([sys.executable,'setup.py','build_ext','--inplace'],cwd=engine,check=True)
 subprocess.run([sys.executable,'native/build_engine.py'],cwd=engine,check=True)
print(engine)
