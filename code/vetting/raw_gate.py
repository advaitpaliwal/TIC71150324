# Event-local linear baselines adapted from hunt4/localize/diff.py event_img.
from pathlib import Path
import sys,json,subprocess,ast
import numpy as np,pandas as pd
from astropy.io import fits
r=Path('work/hunt7');tic=int(sys.argv[1]);x=json.load(open(r/f'res/{tic}.json'));star=pd.read_csv(r/'pool_all.csv').set_index('ID').loc[tic]
P,T0=x['P'],x['T0'];a=(6.674e-11*star.mass*1.989e30*(P*86400)**2/(4*np.pi**2))**(1/3);D=max(.85*P/np.pi*np.arcsin(min(1,star.rad*6.957e8/a)),.03)
ts=f'{tic:016d}';path='/'.join(ts[i:i+4] for i in range(0,16,4));rows=[];missing=[];available=[]
for s in ast.literal_eval(star.sectors):
 for prod in ['spoc','qlp']:
  fn=r/f'lc/{tic}/{prod}_s{s}.fits'
  url=(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{path}/hlsp_tess-spoc_tess_phot_{ts}-s{s:04d}_tess_v1_lc.fits' if prod=='spoc' else f'https://archive.stsci.edu/hlsps/qlp/s{s:04d}/{path}/hlsp_qlp_tess_ffi_s{s:04d}-{ts}_tess_v01_llc.fits')
  if not fn.exists():subprocess.run(['curl','-sf','-m','25','--retry','4','--retry-max-time','120','-o',str(fn),url])
  try:
   with fits.open(fn) as h:
    available.append(f'{prod}_{s}')
    data=h[1].data;t=np.array(data['TIME']);raw=np.array(data['SAP_FLUX']);names=data.names;col='PDCSAP_FLUX' if prod=='spoc' else ('KSPSAP_FLUX' if 'KSPSAP_FLUX' in names else 'DET_FLUX');corr=np.array(data[col]);q=np.array(data['QUALITY'])
   good=np.isfinite(t)&np.isfinite(raw)&np.isfinite(corr)&(q==0);t,raw,corr=t[good],raw[good],corr[good]
   for n in range(int(np.floor((t.min()-T0)/P)),int(np.ceil((t.max()-T0)/P))+1):
    dt=t-(T0+n*P);inside=np.abs(dt)<.4*D;outside=(np.abs(dt)>D)&(np.abs(dt)<.5)
    if inside.sum()<2 or outside.sum()<10 or (dt[outside]<0).sum()<3 or (dt[outside]>0).sum()<3:continue
    result=dict(tic=tic,product=prod,sector=s,event=n,time=T0+n*P,nin=int(inside.sum()))
    for label,flux in [('raw',raw),('corrected',corr)]:
     X=np.vstack([np.ones(outside.sum()),dt[outside]]).T;Y=flux[outside];beta=np.linalg.lstsq(X,Y,rcond=None)[0];res=Y-X@beta;scale=1.4826*np.median(np.abs(res-np.median(res)));ok=np.abs(res)<5*scale+1e-9
     if ok.sum()>8:beta=np.linalg.lstsq(X[ok],Y[ok],rcond=None)[0]
     base=beta[0]+dt[inside]*beta[1];depth=np.mean((base-flux[inside])/base);err=np.std((Y-X@beta)/beta[0])*np.sqrt(1/inside.sum()+1/outside.sum())
     result[label+'_ppm']=float(depth*1e6);result[label+'_err_ppm']=float(err*1e6)
    rows.append(result)
  except Exception as e:missing.append(dict(product=prod,sector=s,error=str(e)[:150]))
d=pd.DataFrame(rows);d.to_csv(r/f'raw_events_{tic}.csv',index=False)
def combine(g,label):
 v=g[label+'_ppm'].values;e=g[label+'_err_ppm'].values;m=np.isfinite(v)&np.isfinite(e)&(e>0);w=1/e[m]**2
 return {'depth_ppm':float(np.sum(v[m]*w)/w.sum()),'err_ppm':float(1/np.sqrt(w.sum())),'events':int(m.sum())} if len(w) else None
out={'tic':tic,'P':P,'T0':T0,'duration_days':D,'baseline_days':.5,'available':available,'missing':missing,'combined':{},'per_sector':{}}
if len(d):
 for prod,g in d.groupby('product'):
  out['combined'][prod]={label:combine(g,label) for label in ['raw','corrected']}
 for (prod,s),g in d.groupby(['product','sector']):out['per_sector'][f'{prod}_{s}']={label:combine(g,label) for label in ['raw','corrected']}
json.dump(out,open(r/f'raw_gate_{tic}.json','w'),indent=2);print(json.dumps(out),flush=True)
