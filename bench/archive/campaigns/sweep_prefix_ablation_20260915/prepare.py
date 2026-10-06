from pathlib import Path
import argparse,os,shutil,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);a=p.parse_args()
package=Path(__file__).resolve().parent;work=Path(os.environ.get('MAPPA_EXPERIMENT_WORK',package/'work')).resolve()
env=dict(os.environ,MAPPA_EXPERIMENT_WORK=str(work))
subprocess.run([sys.executable,str(package.parent/'prefix_competition_20260915/prepare.py'),'--repo',str(a.repo.resolve())],env=env,check=True)
engine=work/'engine';shutil.copy2(engine/'src/mappa/competition.py',engine/'src/mappa/paired_competition.py')
subprocess.run(['patch','-p1','-i',str(package/'paired.patch')],cwd=engine,check=True)
