# all6 fitdata for TIC 71150324: original 3 files (SPOC s5, SPOC s32, QLP s98 DET_FLUX; quality==0) + coverage-lane TESScut SAP s106/107/108 (200 s).
# TESScut mask = sectors/scripts/measure.py: SPOC default bitmask 0b0101001010111111 plus bit 2048 (stray light) -> drops epoch 340 (all 2048, bkg ~1850 vs ~98).
# Epoch 333 has no cadences in the S106 TESScut LC. TESScut TARGET has no error column -> error = per-sector OOT MAD; no crowding correction (raw SAP).
# Detrend: wotan biweight 0.5 d, transit masked 3x dur; keep |phase|<0.15 d; 5/-10 sigma clip. Source group: 0 SPOC1800, 1 SPOC600, 2 QLP200, 3 TESSCUT200.
import numpy as np, json
from astropy.io import fits
from wotan import flatten, transit_mask
HERE='work/deep/71150324/fit/all6'; LC='work/hunt7/lc/71150324'; TC='work/deep/71150324/sectors/tesscut_lc'
R=json.load(open('work/hunt7/res/71150324.json')); P0,T00=R['P'],R['T0']
D0=json.load(open('work/hunt7/raw_gate_71150324.json'))['duration_days']
FILES=[('spoc_s5',f'{LC}/spoc_s5.fits','PDCSAP_FLUX','PDCSAP_FLUX_ERR',0),('spoc_s32',f'{LC}/spoc_s32.fits','PDCSAP_FLUX','PDCSAP_FLUX_ERR',1),('qlp_s98',f'{LC}/qlp_s98.fits','DET_FLUX','DET_FLUX_ERR',2)]
FILES+=[(f'tesscut_s{s}',f'{TC}/s{s}_sap.fits','TARGET',None,3) for s in (106,107,108)]
T,F,E,S,X,G=[],[],[],[],[],[]; log={}
for nm,fn,fc,ec,gid in FILES:
    h=fits.open(fn); d=h[1].data; s=h[0].header['SECTOR']; t0=np.array(d['TIME'],float)
    exp=int(round(h[1].header['TIMEDEL']*86400)) if 'TIMEDEL' in h[1].header else int(round(np.median(np.diff(t0))*86400))
    q=np.array(d['QUALITY']).astype(np.int64)
    good=(q==0) if ec else ((q&(0b0101001010111111|2048))==0)
    noerr=ec is None
    fl0=np.array(d[fc],float); m=good&np.isfinite(fl0)&np.isfinite(t0)&(noerr or np.isfinite(np.array(d[ec],float)))
    t,f=t0[m],fl0[m]; e=np.ones_like(f) if noerr else np.array(d[ec],float)[m]; e=e/np.nanmedian(f); f=f/np.nanmedian(f)
    msk=transit_mask(time=t,period=P0,duration=3*D0,T0=T00)
    tr=flatten(t[~msk],f[~msk],method='biweight',window_length=0.5,return_trend=True)[1]
    gg=np.isfinite(tr); fl=np.interp(t,t[~msk][gg],tr[gg]); f,e=f/fl,e/fl
    r=f-1; sig=1.4826*np.median(np.abs(r[~msk]-np.median(r[~msk])))
    if noerr: e=np.full_like(f,sig)
    ph=((t-T00+0.5*P0)%P0)-0.5*P0; w=(r<5*sig)&(r>-10*sig)&(np.abs(ph)<0.15)
    ep=np.unique(np.round((t[w&(np.abs(ph)<D0/2)]-T00)/P0).astype(int)).tolist()
    T+=[t[w]];F+=[f[w]];E+=[e[w]];S+=[np.full(w.sum(),s)];X+=[np.full(w.sum(),exp)];G+=[np.full(w.sum(),gid)]
    log[nm]=dict(sector=s,group=gid,exp_s=exp,kept=int(w.sum()),in_transit_pts=int(((np.abs(ph)<D0/2)&w).sum()),epochs_with_intransit=ep,oot_mad_ppm=round(sig*1e6),med_err_ppm=round(float(np.median(e))*1e6),
        quality_mask=('QUALITY==0' if ec else 'QUALITY & (0b0101001010111111|2048)==0 (measure.py)'),err_column=('none: OOT MAD' if noerr else ec),flux_column=fc,n_quality_rejected=int((~good).sum()))
    print(nm,log[nm],flush=True)
np.savez(HERE+'/fitdata.npz',t=np.concatenate(T),f=np.concatenate(F),e=np.concatenate(E),s=np.concatenate(S),x=np.concatenate(X),g=np.concatenate(G))
json.dump(dict(P0=P0,T0_0=T00,dur0_d=D0,files=log),open(HERE+'/prep_log.json','w'),indent=1)
