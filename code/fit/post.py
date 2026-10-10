# Post-process chain -> result json + folded model-over-data PNG (TIC 71150324). usage: post.py [chain.npz] [tag]
import numpy as np, sys, json, batman
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
CH = sys.argv[1] if len(sys.argv) > 1 else 'chain.npz'; TAG = sys.argv[2] if len(sys.argv) > 2 else ''
sys.argv = ['x']; sys.path.insert(0, 'work/deep/71150324/fit/all6'); import fit as Fm
HERE = Fm.HERE; np.random.seed(1)
c = np.load(HERE+'/'+CH, allow_pickle=True); ch = c['chain']; N = len(ch)
T0, P, k, b, rho, q1, q2, f0 = ch.T[:8]
a = Fm.a_rs(rho, P); inc = np.degrees(np.arccos(b/a)); u1 = 2*np.sqrt(q1)*q2; u2 = np.sqrt(q1)*(1-2*q2)
T14 = P/np.pi*np.arcsin(np.sqrt(np.clip((1+k)**2-b**2, 0, None))/(a*np.sin(np.radians(inc))))*24
def bp(i):
    p = batman.TransitParams(); p.t0=0; p.per=P[i]; p.rp=k[i]; p.a=a[i]; p.inc=inc[i]; p.ecc=0; p.w=90; p.limb_dark='quadratic'; p.u=[u1[i], u2[i]]; return p
sel = np.random.choice(N, min(N, 4000), replace=False)
dep = np.array([1-batman.TransitModel(bp(i), np.array([0.])).light_curve(bp(i))[0] for i in sel])*1e6
st = Fm.st
Rs = np.random.normal(float(st['srad']), float(st['srad_e']), N); Ms = np.random.normal(float(st['mass']), float(st['mass_e']), N)
Teff = np.random.normal(float(st['teff']), float(st['teff_e']), N); L = np.random.normal(float(st['lum']), float(st['lum_e']), N)
Rp = k*Rs*109.076
sma = (Ms*(P/365.25)**2)**(1/3)
Teq = Teff*np.sqrt(1/(2*a))     # Bond albedo A_B=0, full day-night redistribution (f=1/4); a/Rs from chain
S = L/sma**2                    # S_earth: TIC luminosity / (Kepler III a with TIC mass)^2
def q(x): m, lo, hi = np.percentile(x, [50, 15.865, 84.135]); return [float(m), float(hi-m), float(m-lo)]
out = {'quantiles_format': '[median, +1sigma, -1sigma] (16/50/84 pct)'}
for nm, x in [('P_d', P), ('T0_BJD_TDB', T0+2457000.), ('depth_ppm', dep), ('T14_h', T14), ('k_RpRs', k), ('a_Rs', a), ('b', b), ('inc_deg', inc),
              ('rho_star_gcc', rho), ('u1', u1), ('u2', u2), ('Rp_Re', Rp), ('Teq_K_AB0_f1over4', Teq), ('S_Se', S), ('a_au', sma)]:
    out[nm] = q(x)
for i, nm in enumerate(c['names'][8:]): out['jitter_ppm_'+nm[4:]] = q(np.exp(ch[:, 8+i])*1e6)
tau = c['tau']; ns = int(c['nsteps']); nw = int(c['nw']); burn = int(c['burn'])
out['convergence'] = dict(nsteps=ns, nwalkers=nw, burn=burn, tau=dict(zip(c['names'].tolist(), np.round(tau, 1).tolist())),
    nsteps_over_max_tau=float(ns/np.nanmax(tau)), postburn_steps_over_max_tau=float((ns-burn)/np.nanmax(tau)), ess_min_postburn=float((ns-burn)*nw/np.nanmax(tau)), ess_post_burn=dict(zip(c['names'].tolist(), np.round((ns-burn)*nw/tau).tolist())),
    rhat_note='Gelman-Rubin computed treating emcee walkers as chains: heuristic only, walkers are correlated (not independent chains)', rhat_walkers_heuristic=dict(zip(c['names'].tolist(), np.round(c['rhat'], 4).tolist())), acc_mean=float(np.mean(c['acc'])), N_thinned=N)
out['b_95pct_upper'] = float(np.percentile(b, 95)); out['corr_T0_P'] = float(np.corrcoef(T0, P)[0, 1]); out['epoch_ref_n_from_hunt7_T0'] = int(c['TREF_N'])
out['ld_prior'] = [float(c['U1']), float(c['U2']), 0.1]; out['ld_caveat'] = 'Claret 2017 table25 ATLAS LSM, Z=0 xi=2, bilinear in Teff (3500-3750 K) and logg (4.5-5.0) at TIC Teff 3586 K logg 4.744; Gaussian prior sigma 0.1 on u1,u2'; out['rho_prior'] = [Fm.RHO, Fm.E_RHO] if Fm.RHO_PRIOR else 'uniform 0.5-40 g/cc (no TIC prior)'
out['rho_caveat'] = 'rho_star posterior is conditioned on the TIC density prior; agreement with TIC is not independent validation'; out['teq_assumption'] = 'Bond albedo 0, full redistribution f=1/4, Teq=Teff*sqrt(Rs/2a)'
json.dump(out, open(HERE+'/result%s.json' % TAG, 'w'), indent=1)
for kx, v in out.items(): print(kx, v)
if TAG: sys.exit()

med = np.median(ch, 0); Pm, Tm = np.median(P), np.median(T0)
ph = ((Fm.T - Tm + 0.5*Pm) % Pm - 0.5*Pm)*24
fn = Fm.F/np.median(f0); res = Fm.F - Fm.model(med)
fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
for A in ax:
    for s in ('top', 'right'): A.spines[s].set_visible(False)
    A.grid(axis='y', c='0.92', lw=0.6)
tt = np.linspace(-0.12, 0.12, 3000)
pm = bp(0); pm.per=Pm; pm.rp=np.median(k); pm.a=np.median(a); pm.inc=np.median(inc); pm.u=[np.median(u1), np.median(u2)]
for i in np.random.choice(N, 50):
    ax[0].plot(tt*24, (batman.TransitModel(bp(i), tt).light_curve(bp(i))-1)*1e6, c='#2a78d6', alpha=0.06, lw=1)
w = 2.5
bins = np.arange(-w, w+1e-9, 0.1); bc = (bins[1:]+bins[:-1])/2; idx = np.digitize(ph, bins)
def bstat(y):
    m = np.array([np.mean(y[idx == j]) if (idx == j).sum() > 3 else np.nan for j in range(1, len(bins))])
    e = np.array([np.std(y[idx == j])/np.sqrt((idx == j).sum()) if (idx == j).sum() > 3 else np.nan for j in range(1, len(bins))]); return m, e
mb, eb = bstat(fn); rb, re_ = bstat(res)
ax[0].plot(ph, (fn-1)*1e6, '.', c='0.82', ms=1.5, zorder=0, label='TESS S5/S32 SPOC, S98 QLP, S106-108 TESScut SAP (unbinned)')
ax[0].errorbar(bc, (mb-1)*1e6, eb*1e6, fmt='o', c='0.15', ms=3.5, lw=1, label='data, 6-min phase bins')
# Draw the same finite exposures as the fit; pale draws above are instantaneous.
for ex, ss, color, style in zip(Fm.EXPS[:3], Fm.SS[:3], ['#2a78d6', '#eb6834', '#1baf7a'], ['-', '--', ':']):
    integrated = batman.TransitModel(pm, tt, supersample_factor=ss, exp_time=ex/86400.).light_curve(pm)
    ax[0].plot(tt*24, (integrated-1)*1e6, c=color, ls=style, lw=1.6, label='model integrated over %d s' % ex)
model_binned, _ = bstat(Fm.model(med)/np.median(f0))
ax[0].plot(bc, (model_binned-1)*1e6, 'x', c='#52514e', ms=4, label='median model, binned like data')
ax[0].set_ylim(-3.0*np.median(dep), 1.5*np.median(dep)); ax[0].set_ylabel('relative flux [ppm]'); ax[0].legend(fontsize=7, loc='lower left', frameon=False)
ax[0].set_title('TIC 71150324 all6 (S5,32,98,106,107,108)   P = %.6f d   depth %.0f ppm   Rp = %.2f R$_\\oplus$' % (Pm, np.median(dep), np.median(Rp)), fontsize=10)
ax[1].errorbar(bc, rb*1e6, re_*1e6, fmt='o', c='0.15', ms=3, lw=1); ax[1].axhline(0, c='#2a78d6', lw=1.6)
ax[1].set_ylabel('resid. [ppm]'); ax[1].set_xlabel('hours from mid-transit'); ax[1].set_xlim(-w, w)
plt.tight_layout(); plt.savefig(HERE+'/folded_model_TIC71150324_all6.png', dpi=130)
