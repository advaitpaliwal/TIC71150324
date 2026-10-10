# usage: run1.py TIC R M "s1,s2,..."   -> downloads TESS-SPOC + QLP FFI LCs, one product per sector, TLS
import sys, os, json, warnings, subprocess, numpy as np
warnings.filterwarnings('ignore')
import lightkurve as lk
from transitleastsquares import transitleastsquares
tic=int(sys.argv[1]); R=float(sys.argv[2]); M=float(sys.argv[3]); secs=[int(s) for s in sys.argv[4].split(',')]
out=f'res/{tic}.json'
if os.path.exists(out) or str(tic) in open('skip.txt').read().split(): sys.exit()
ts=f'{tic:016d}'; p=f'{ts[0:4]}/{ts[4:8]}/{ts[8:12]}/{ts[12:16]}'
def get(url,fn):
    if not os.path.exists(fn):
        r=subprocess.run(['curl','-sf','-m','60','--retry','4','--retry-max-time','180','-o',fn,url])
        if r.returncode!=0:
            if os.path.exists(fn): os.remove(fn)
            return None
    return fn
os.makedirs(f'lc/{tic}',exist_ok=True)
parts=[]; used=[]
try:
    for s in secs:
        lc=None
        fn=get(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{p}/hlsp_tess-spoc_tess_phot_{ts}-s{s:04d}_tess_v1_lc.fits',f'lc/{tic}/spoc_s{s}.fits')
        if fn:
            try: lc=lk.read(fn,flux_column='pdcsap_flux',quality_bitmask='hard'); src='TESS-SPOC'
            except Exception: lc=None
        if lc is None:
            fn=get(f'https://archive.stsci.edu/hlsps/qlp/s{s:04d}/{p}/hlsp_qlp_tess_ffi_s{s:04d}-{ts}_tess_v01_llc.fits',f'lc/{tic}/qlp_s{s}.fits')
            if fn:
                for col in ('kspsap_flux','det_flux'):
                    try: lc=lk.read(fn,flux_column=col,quality_bitmask='hard'); src='QLP'; break
                    except Exception: lc=None
        if lc is None: continue
        lc=lc.remove_nans().normalize()
        if len(lc)<100: continue
        cad=np.nanmedian(np.diff(lc.time.value))
        wl=int(0.75/cad)//2*2+1
        lc=lc.flatten(window_length=wl).remove_outliers(sigma_upper=3,sigma_lower=20)
        lc=lc.bin(time_bin_size=max(cad,30/1440)).remove_nans()
        parts.append(lc); used.append((s,src,round(cad*1440,1)))
    if len(parts)<(2 if len(secs)==2 else 3): raise RuntimeError(f'only {len(parts)} sectors')
    t=np.concatenate([x.time.value for x in parts]); f=np.concatenate([np.asarray(x.flux.value,float) for x in parts])
    sid=np.concatenate([np.full(len(x),u[0]) for x,u in zip(parts,used)])
    o=np.argsort(t); t,f,sid=t[o],f[o],sid[o]
    np.savez(f'lc/{tic}.npz',t=t,f=f,s=sid)
    m=transitleastsquares(t,f)
    r=m.power(period_min=1.0,period_max=min(12,(t.max()-t.min())/2),R_star=R,M_star=M,R_star_min=0.1,R_star_max=1.2,M_star_min=0.1,M_star_max=1.2,oversampling_factor=1,duration_grid_step=1.1,use_threads=1,show_progress_bar=False)
    d=dict(tic=tic,used=used,SDE=float(r.SDE),P=float(r.period),T0=float(r.T0),depth=float(1-r.depth),dur=float(r.duration),snr=float(r.snr),oe=float(r.odd_even_mismatch),ntr=int(r.distinct_transit_count),rp_rs=float(r.rp_rs),std=float(np.std(f)))
except Exception as e:
    d=dict(tic=tic,used=used,err=repr(e)[:300])
json.dump(d,open(out,'w'))
print(d)
