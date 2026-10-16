# Rev1 FPP for TIC 71150324 (AAS82162). Copied from prerev/71150324/fpp/run.py (itself a copy of the submitted joint run).
# Frozen: TRICERATOPS 1.1.0 venv, pickle, Gaia background, MOLUSC file, filt, N=1e6, exptime=30/1440, nsamples=20, seeds 1000+k, 10-min phase bins, vet duration.
# Changed for Revision 1: LC = all six sectors (deep/71150324/fit/all6/fitdata.npz: S5, S32, S98, S106-108, wotan biweight with transit mask, |phase|<0.15 d),
# ephemeris = published fit (P, T0 from fit/all6/result.json), tdepth = fitted depth 2277.375 ppm.
# Modes: full = no neighbours removed; sel = remove TICs in drop_sel.txt.
import sys, json, pickle, os, time, numpy as np, warnings; warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); PH = '~/Scratch/planet-hunt/'; D = PH + 'deep/71150324/fpp/'
sys.path.insert(0, D + 'scripts'); sys.path.insert(0, PH + 'hunt5b/trc/vendor')
import gaia_patch
import triceratops.triceratops as tr
import triceratops.marginal_likelihoods as ml, types
class _F:
    def __init__(s, d): s.d = d
    def __getitem__(s, k):
        if isinstance(k, str): return types.SimpleNamespace(values=s.d[k].to_numpy(copy=True))
        return _F(s.d[np.asarray(k)])
_rc = ml.read_csv
ml.read_csv = lambda f, *a, **kw: _F(_rc(f, *a, **kw)) if 'molusc' in str(f) else _rc(f, *a, **kw)
mode, k = sys.argv[1], int(sys.argv[2]); assert mode in ('full', 'sel')
R = json.load(open(PH + 'deep/71150324/fit/all6/result.json'))
P = R['P_d'][0]; T0 = R['T0_BJD_TDB'][0] - 2457000; depth = R['depth_ppm'][0] * 1e-6
z = np.load(PH + 'deep/71150324/fit/all6/fitdata.npz'); t, f = z['t'], z['f']
tg = pickle.load(open(D + 'out/target_71150324.pkl', 'rb'))
drop = [] if mode == 'full' else open(HERE + '/drop_sel.txt').read().split()
before = len(tg.stars)
if drop:
    keep = ~tg.stars['ID'].astype(str).isin(drop).values
    tg.pix_coords = [pc[keep] for pc in tg.pix_coords]; tg.remove_star(drop); tg.stars = tg.stars.reset_index(drop=True)
if tg.trilegal_fname: tg.trilegal_fname = D + '71150324_gaia_background.csv'
ph = t - T0 - np.round((t - T0) / P) * P; w = np.abs(ph) < 0.15
edges = np.arange(-0.15, 0.15 + 1e-9, 10 / 1440); idx = np.digitize(ph[w], edges)
tb = np.array([ph[w][idx == i].mean() for i in np.unique(idx)]); fb = np.array([f[w][idx == i].mean() for i in np.unique(idx)])
dur = json.load(open(PH + 'hunt7/vet_71150324.json'))['dur_h'] / 24
oot = np.abs(tb) > dur; err = float(np.std(fb[oot])); dbin = float(np.median(fb[oot]) - fb[np.abs(tb) < dur / 2 * 0.8].mean())
tg.calc_depths(tdepth=depth, all_ap_pixels=tg._aps)
np.random.seed(1000 + k); t0 = time.time()
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=err, P_orb=P, contrast_curve_file=None, filt='TESS', N=1000000, parallel=True,
              exptime=30 / 1440, nsamples=20, verbose=0, molusc_file=D + 'molusc_out/molusc_for_trc.csv')
tag = f'{mode}_{k}'; tg.probs.to_csv(HERE + f'/out/probs_{tag}.csv', index=False)
json.dump(dict(mode=mode, run=k, seed=1000 + k, dropped=drop, n_stars=len(tg.stars), n_before=before, P=P, T0_btjd=T0, depth_ppm=depth * 1e6,
               binned_box_depth_ppm=dbin * 1e6, err_ppm=err * 1e6, nbins=len(tb), duration_days=dur, sectors=sorted(set(int(s) for s in z['s'])),
               FPP=float(tg.FPP), NFPP=float(tg.NFPP), sec=time.time() - t0), open(HERE + f'/out/res_{tag}.json', 'w'))
print(tag, tg.FPP, tg.NFPP, round(time.time() - t0), flush=True)
