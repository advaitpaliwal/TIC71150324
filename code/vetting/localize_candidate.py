# Reuse hunt5a difference imaging; only substitute target and ephemeris.
from pathlib import Path
import sys,json,subprocess,os,time
import pandas as pd
r=Path(__file__).resolve().parent;tic=int(sys.argv[1]);py=str(r.parent/'hunt4/.venv/bin/python')
x=json.loads((r/f'raw_gate_{tic}.json').read_text());star=pd.read_csv(r/'pool_all.csv').set_index('ID').loc[tic]
root=r/f'localize/{tic}';root.mkdir(exist_ok=True);(root/'out').mkdir(exist_ok=True)
src=r.parent/'hunt5a/localize/22205102'
s=(src/'common.py').read_text().replace('cut/s{s}_*/tess-s','cut/s{s}_*/tesscut/tess-s');lines=s.splitlines();lines[4]=f"P,T0,DUR={x['P']},{x['T0']},{x['duration_days']}";lines[5]=f'RA0,DE0={star.ra},{star.dec}';(root/'common.py').write_text('\n'.join(lines)+'\n')
(root/'diff.py').write_text((src/'diff.py').read_text().replace('g[g.sep<=60]', 'g[g.sep<=120]'));(root/'dl.py').write_text((src/'dl.py').read_text().replace('TIC 22205102',f'TIC {tic}'))
end=float(os.environ.get('HUNT_DEADLINE',json.loads((r/'pool_counts.json').read_text())['start_unix']+87*60))
sectors=sorted({int(k.split('_')[1]) for k in x['per_sector']})
env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','MPLCONFIGDIR':str(r/'mplcache'),'XDG_CACHE_HOME':str(r/'cache'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
for sector in sectors:
 if time.time()>end:break
 cut=root/f'cut/s{sector}_stamp';cut.mkdir(parents=True,exist_ok=True)
 with open(root/f'dl_s{sector}.log','w') as f:
  try:p=subprocess.run([py,'-I',str(root/'dl.py'),str(sector),str(cut)],cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=min(240,max(1,end-time.time())))
  except subprocess.TimeoutExpired:print('DOWNLOAD_TIMEOUT',tic,sector,flush=True);continue
 if p.returncode:print('DOWNLOAD_ERROR',tic,sector,p.returncode,flush=True);continue
 with open(root/f'diff_s{sector}.log','w') as f:
  try:p=subprocess.run([py,'-I',str(root/'diff.py'),str(sector)],cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=min(300,max(1,end-time.time())))
  except subprocess.TimeoutExpired:print('LOCALIZE_TIMEOUT',tic,sector,flush=True);continue
 print('LOCALIZED',tic,sector,p.returncode,flush=True)
 result=root/f'out/diff_s{sector}.json'
 if p.returncode==0 and result.exists():
  im=json.loads(result.read_text());target=min(im['stars'],key=lambda v:v['sep'])
  neighbours=[v for v in im['stars'] if v['source']!=target['source']]
  near=[v for v in neighbours if v['sep']<=60]
  bad=[v for v in neighbours if v['chi2']-target['chi2']<=-25]
  if bad or (not near and (im['best']['sep_from_target_arcsec']>21 or any(v['chi2']<target['chi2'] for v in neighbours))):
   print('SECTOR_GATE_C_FAIL',tic,sector,flush=True);break
print('LOCALIZATION_FINISHED',tic,flush=True)
