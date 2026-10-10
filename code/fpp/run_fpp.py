# Adapted from ../../triceratops/{setup_target,run_fpp}.py for FFI-only stars.
# usage: run_fpp.py TIC P T0 RUN_INDEX
import sys, json, pickle, os, time, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
os.chdir('work/deep/71150324/fpp')
sys.path.insert(0,'work/deep/71150324/fpp/scripts')
sys.path.insert(0,'work/hunt5b/trc/vendor')
import gaia_patch
import triceratops.triceratops as tr
import triceratops.marginal_likelihoods as ml, types
class _F:
    def __init__(s,d): s.d=d
    def __getitem__(s,k):
        if isinstance(k,str): return types.SimpleNamespace(values=s.d[k].to_numpy(copy=True))
        return _F(s.d[np.asarray(k)])
_rc=ml.read_csv
ml.read_csv=lambda f,*a,**kw: _F(_rc(f,*a,**kw)) if 'molusc' in str(f) else _rc(f,*a,**kw)
tic, P, T0, k = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
z = np.load('work/hunt7/lc/71150324.npz'); t, f, sid = z['t'], z['f'], z['s']
secs = np.unique(sid).astype(int)
pk = f'out/target_{tic}.pkl'
if not os.path.exists(pk):
    print('building target',flush=True); tg = tr.target(ID=tic, sectors=secs); print('target ok',flush=True)
    aps = []
    for i, s in enumerate(secs):
        fn = f'work/hunt7/lc/{tic}/spoc_s{s}.fits'
        ap = None
        if os.path.exists(fn):
            with fits.open(fn) as h:
                a = h['APERTURE'].data.astype(int); c0 = h['APERTURE'].header['CRVAL1P']; r0 = h['APERTURE'].header['CRVAL2P']
            rows, cols = np.nonzero(np.bitwise_and(a, 2))
            if len(rows): ap = np.column_stack([cols + c0, rows + r0])
        if ap is None:  # QLP: no aperture in file; 5x5 centred on target (TRICERATOPS default)
            tp = np.round(tg.pix_coords[i][0])
            ap = np.array([np.repeat(np.arange(tp[0]-2, tp[0]+3), 5), np.tile(np.arange(tp[1]-2, tp[1]+3), 5)]).T
        aps.append(ap)
    tg._aps = aps
    pickle.dump(tg, open(pk, 'wb'))
tg = pickle.load(open(pk, 'rb'))
# Same removal pattern as deep/375939912/fpp/scripts/run_fpp.py.
drop = pd.read_csv('exclusions.csv', dtype={'tic': str}).query('excluded')['tic'].tolist()
keep = ~tg.stars['ID'].astype(str).isin(drop).values
assert len(drop) == 7 and (~keep).sum() == 7 and str(tic) not in drop
before = len(tg.stars)
tg.pix_coords = [pc[keep] for pc in tg.pix_coords]
tg.remove_star(drop)
tg.stars = tg.stars.reset_index(drop=True)
assert len(tg.stars) == before - 7
assert all(len(pc) == len(tg.stars) for pc in tg.pix_coords)
assert not tg.stars['ID'].astype(str).isin(drop).any()
# Point any cached background filename at the byte-identical local copy.
if tg.trilegal_fname:
    tg.trilegal_fname = os.path.abspath('71150324_gaia_background.csv')
tg.stars.to_csv(f'out/stars_after_{k}.csv', index=False)

ph = t - T0 - np.round((t - T0) / P) * P
w = np.abs(ph) < 0.25
edges = np.arange(-0.25, 0.25 + 1e-9, 10 / 1440)
idx = np.digitize(ph[w], edges)
tb = np.array([ph[w][idx == i].mean() for i in np.unique(idx)])
fb = np.array([f[w][idx == i].mean() for i in np.unique(idx)])
dur = float(sys.argv[5]) if len(sys.argv) > 5 else json.load(open('work/hunt7/vet_71150324.json'))['dur_h'] / 24
oot = np.abs(tb) > dur
err = float(np.std(fb[oot]))
depth = float(np.median(fb[oot]) - fb[np.abs(tb) < dur / 2 * 0.8].mean())
original = json.load(open('work/hunt7/fpp/71150324/out/res_71150324_0.json'))
assert depth * 1e6 == original['depth_ppm'] and err * 1e6 == original['err_ppm']
tg.calc_depths(tdepth=depth, all_ap_pixels=tg._aps)
np.random.seed(1000 + k); t0 = time.time()
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=err, P_orb=P, contrast_curve_file=None, filt='TESS',
              N=1000000, parallel=True, exptime=30 / 1440, nsamples=20, verbose=0, molusc_file="molusc_out/molusc_for_trc.csv")
tag = f'{tic}_{k}'
tg.probs.to_csv(f'out/probs_{tag}.csv', index=False)
json.dump(dict(tic=tic, run=k, P=P, T0=T0, depth_ppm=depth * 1e6, err_ppm=err * 1e6, nbins=len(tb),
               duration_days=dur, seed=1000+k, N=1000000, dropped=drop, stars_before=before, stars_after=len(tg.stars), pixel_shapes=[list(pc.shape) for pc in tg.pix_coords], molusc=os.path.abspath('molusc_out/molusc_for_trc.csv'), FPP=float(tg.FPP), NFPP=float(tg.NFPP), sec=time.time() - t0), open(f'out/res_{tag}.json', 'w'))
print(tag, tg.FPP, tg.NFPP, round(time.time() - t0))
