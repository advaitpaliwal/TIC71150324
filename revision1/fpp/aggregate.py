# Aggregate rev1 runs: FPP/NFPP stats, scenario means, per-star NFPP shares joined to six-sector pixel tests.
import json, glob, sys, numpy as np, pandas as pd
H = '~/Scratch/planet-hunt/'; O = H + 'revision-82162/fpp/'
loc = pd.read_csv(H + 'prerev/71150324/fpp/exclusions_prerev.csv', dtype={'tic': str, 'gaia': 'Int64'})
ps = pd.read_csv(H + 'prerev/71150324/loc/dchi2_sum_6sectors_persector.csv')
nfpp_sc = ['NTP', 'NEB', 'NEBx2P']
for mode in ('full', 'sel'):
    R = [json.load(open(f)) for f in sorted(glob.glob(O + f'out/res_{mode}_*.json'))]
    if not R: continue
    F = np.array([r['FPP'] for r in R]); N = np.array([r['NFPP'] for r in R])
    pr = [pd.read_csv(O + f"out/probs_{mode}_{r['run']}.csv", dtype={'ID': str}) for r in R]
    sc = pd.concat([p.groupby('scenario').prob.sum() for p in pr], axis=1).mean(axis=1)
    shares = pd.concat([p[p.scenario.isin(nfpp_sc)].groupby('ID').prob.sum() / p[p.scenario.isin(nfpp_sc)].prob.sum() for p in pr], axis=1)
    sh = pd.DataFrame(dict(share_mean=shares.mean(axis=1), share_min=shares.min(axis=1), share_max=shares.max(axis=1))).sort_values('share_mean', ascending=False)
    sh = sh[sh.share_mean > 0].reset_index().rename(columns={'ID': 'tic'}).merge(loc, on='tic', how='left')
    gaia = dict(zip(loc.tic, loc.gaia))
    def rule(row):
        if pd.isna(row.gaia): return 'no pixel test'
        d = ps[ps.source == int(row.gaia)].delta_chi2
        return 'yes' if (d.sum() >= 25 and (d > -25).all()) else 'no'
    sh['rule_met'] = sh.apply(rule, axis=1)
    sh.to_csv(O + f'perstar_nfpp_{mode}.csv', index=False)
    out = dict(mode=mode, n=len(R), runs=[dict(run=r['run'], seed=r['seed'], FPP=r['FPP'], NFPP=r['NFPP']) for r in R],
               FPP_mean=F.mean(), FPP_sd=F.std(ddof=1) if len(F) > 1 else None, FPP_min=F.min(), FPP_max=F.max(),
               n_le_0p035=int((F <= 0.035).sum()), n_lt_0p015=int((F < 0.015).sum()),
               NFPP_mean=N.mean(), NFPP_sd=N.std(ddof=1) if len(N) > 1 else None,
               n_stars=R[0]['n_stars'], dropped=R[0]['dropped'], depth_ppm=R[0]['depth_ppm'], binned_box_depth_ppm=R[0]['binned_box_depth_ppm'],
               err_ppm=R[0]['err_ppm'], nbins=R[0]['nbins'], sectors=R[0]['sectors'], P=R[0]['P'], T0_btjd=R[0]['T0_btjd'],
               scenarios=sc.to_dict())
    json.dump(out, open(O + f'summary_{mode}.json', 'w'), indent=1, default=float)
    print(mode, len(R), 'FPP %.4f +- %s [%.4f, %.4f] le0.035 %d  NFPP %.4g' % (F.mean(), out['FPP_sd'], F.min(), F.max(), out['n_le_0p035'], N.mean()))
    print(sh.to_string()); print(sc.sort_values(ascending=False).head(10).to_string())
