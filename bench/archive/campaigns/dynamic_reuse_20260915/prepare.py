from pathlib import Path
import argparse,os,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--repo',required=True,type=Path);p.add_argument('--reference',action='store_true');a=p.parse_args()
package=Path(__file__).resolve().parent
work=Path(os.environ.get('MAPPA_EXPERIMENT_WORK',package/'work')).resolve()
env=dict(os.environ,MAPPA_EXPERIMENT_WORK=str(work))
subprocess.run([sys.executable,str(package.parent/'prefix_competition_20260915/prepare.py'),'--repo',str(a.repo.resolve())],env=env,check=True)
if not a.reference:
 subprocess.run(['patch','-p1','-i',str(package/'optimization.patch')],cwd=work/'engine',check=True)
