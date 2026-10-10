# Glue (adapted from deep/71150324/fit/fit.py): all6 batman + emcee fit for TIC 71150324 (adds TESScut SAP S106-108).
# Shared T0,P,k,b,rho,LD,f0; separate white-noise jitter per source group (SPOC1800, SPOC600, QLP200, TESSCUT200). No dilution term.
# Input: fitdata.npz (prep.py): t BTJD (TDB), f detrended, e, s sector, x native exposure s (1800/600/200 from TIMEDEL).
import numpy as np, emcee, batman, json, sys, time
from scipy.optimize import minimize
HERE = 'work/deep/71150324/fit/all6'
np.random.seed(42)
import os
RHO_PRIOR = os.environ.get('NO_RHO_PRIOR') != '1'   # NO_RHO_PRIOR=1 -> uniform rho, transit-only density
OUT = 'chain.npz' if RHO_PRIOR else 'chain_norho.npz'
R = json.load(open('work/hunt7/res/71150324.json')); P0, T00 = R['P'], R['T0']
st = [s for s in json.load(open('work/deep/71150324/fit/_sources/exofop/tic_exofop.json'))['stellar_parameters'] if s.get('snotes') == 'TIC v8.2'][0]
RHO, E_RHO = float(st['dens']), float(st['dens_e']); TEFF, LOGG = float(st['teff']), float(st['logg'])
# Claret 2017 TESS quadratic LD (table25, ATLAS LSM, Z=0, xi=2), bilinear in Teff (3500-3750) and logg (4.5-5.0) at TIC Teff/logg.
g = {}
for ln in open('work/deep/71150324/fit/_sources/claret/claret_raw.tsv'):
    c = ln.split('\t')
    try: lg, T_, Z, xi, a_, b_ = map(float, c[1:7])
    except Exception: continue
    if Z == 0 and xi == 2: g[(T_, lg)] = (a_, b_)
y_ = (LOGG-4.5)/0.5; x_ = (TEFF-3500.)/250.
def bil(j): return (1-x_)*((1-y_)*g[(3500., 4.5)][j] + y_*g[(3500., 5.0)][j]) + x_*((1-y_)*g[(3750., 4.5)][j] + y_*g[(3750., 5.0)][j])
U1 = bil(0); U2 = bil(1); E_U = 0.1

D = np.load(HERE+'/fitdata.npz'); T, F, E, S, X, GID = D['t'], D['f'], D['e'], D['s'], D['x'], D['g']
GLAB = ['spoc1800', 'spoc600', 'qlp200', 'tesscut200']; NG = 4; grp = GID.astype(int)
EXPS = [int(np.round(X[grp == i][0])) for i in range(NG)]
TREF_N = int(np.round((np.median(T) - T00)/P0)); TREF0 = T00 + TREF_N*P0
p = batman.TransitParams(); p.t0=TREF0; p.per=P0; p.rp=0.025; p.a=10.; p.inc=89.; p.ecc=0.; p.w=90.; p.limb_dark='quadratic'; p.u=[U1, U2]
SS = [max(1, int(np.ceil(EXPS[i]/120))) for i in range(NG)]
MODS = [batman.TransitModel(p, T[grp == i], supersample_factor=SS[i], exp_time=EXPS[i]/86400.) for i in range(NG)]
def a_rs(rho, P): G = 6.674e-8; return ((G*rho*(P*86400)**2)/(3*np.pi))**(1/3)
names = ['T0', 'P', 'k', 'b', 'rho', 'q1', 'q2', 'f0'] + ['lns_' + l for l in GLAB]
def unpack(th):
    T0, P, k, b, rho, q1, q2, f0 = th[:8]
    a = a_rs(rho, P); inc = np.degrees(np.arccos(min(b/a, 1)))
    return T0, P, k, b, rho, a, inc, 2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2), f0
def model(th):
    T0, P, k, b, rho, a, inc, u1, u2, f0 = unpack(th)
    p.t0=T0; p.per=P; p.rp=k; p.a=a; p.inc=inc; p.u=[u1, u2]
    out = np.empty_like(F)
    for i in range(NG): out[grp == i] = MODS[i].light_curve(p)
    return out*f0
def lnprob(th):
    T0, P, k, b, rho, q1, q2, f0 = th[:8]; ls = np.asarray(th[8:])
    if not (0 < k < 0.3 and 0 <= b < 1+k and 0.5 < rho < 40 and 0 < q1 < 1 and 0 < q2 < 1 and np.all((ls > -14) & (ls < -4))): return -np.inf
    if abs(T0 - TREF0) > 0.05 or abs(P - P0) > 0.001: return -np.inf
    _, _, _, _, _, a, inc, u1, u2, _ = unpack(th)
    if b >= a: return -np.inf
    lp = (-0.5*((rho-RHO)/E_RHO)**2 if RHO_PRIOR else 0.) - 0.5*((u1-U1)/E_U)**2 - 0.5*((u2-U2)/E_U)**2
    var = E**2 + np.exp(2*ls)[grp]; r = F - model(th)
    return lp - 0.5*np.sum(r**2/var + np.log(var))

if __name__ == '__main__':
    print('LD', U1, U2, 'rho', RHO, E_RHO, 'exps', EXPS, 'supersample', SS, 'N', len(T), 'N per exp', [int((grp == i).sum()) for i in range(NG)], 'TREF_N', TREF_N, 'TREF0', TREF0, flush=True)
    q1_0 = (U1+U2)**2; q2_0 = U1/(2*(U1+U2))
    best = np.array([TREF0, P0, np.sqrt(R['depth']), 0.3, RHO, q1_0, q2_0, 1.0] + [-9]*NG)
    t_ = time.time(); [lnprob(best) for _ in range(50)]; print('sec per lnprob', (time.time()-t_)/50, flush=True)
    nll = lambda x: -lnprob(x)
    for _ in range(4):
        r = minimize(nll, best, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-9, fatol=1e-3)); best = r.x
    print('MAP', dict(zip(names, best)), -r.fun, flush=True)
    nd = len(best); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 3*nd
    sc = np.array([5e-4, 2e-7, 1e-3, 0.05, 0.2, 0.02, 0.02, 1e-5] + [0.05]*NG)
    pos = []
    while len(pos) < nw:
        c = best + sc*np.random.randn(nd)
        if np.isfinite(lnprob(c)): pos.append(c)
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 15000
    sam = emcee.EnsembleSampler(nw, nd, lnprob)
    for i, _ in enumerate(sam.sample(np.array(pos), iterations=nsteps)):
        if (i+1) % 2000 == 0: print('step', i+1, 'tau', np.round(sam.get_autocorr_time(tol=0), 0).tolist(), flush=True)
    tau = sam.get_autocorr_time(tol=0); print('tau', dict(zip(names, np.round(tau, 1))))
    burn = nsteps//3
    ch = sam.get_chain(discard=burn, thin=10, flat=True)
    pc = sam.get_chain()[burn:]; W_ = pc.var(0, ddof=1).mean(0); B_ = pc.mean(0).var(0, ddof=1)*pc.shape[0]
    n_ = pc.shape[0]; rhat = np.sqrt(((n_-1)/n_*W_ + B_/n_)/W_)
    np.savez(HERE+'/'+OUT, chain=ch, tau=tau, names=names, nsteps=nsteps, burn=burn, nw=nw, rhat=rhat,
             acc=sam.acceptance_fraction, TREF_N=TREF_N, U1=U1, U2=U2, best=best, lnp=sam.get_log_prob(discard=burn, thin=10, flat=True))
    print('rhat', dict(zip(names, np.round(rhat, 4))), 'acc', np.mean(sam.acceptance_fraction), 'n', len(ch))
