# The gap from a hitting time (after Roberts and Shortland's hazard rate view of boundary hitting). Brownian motion on SU(2) = S^3 with
# generator (1/2) Laplacian, started at the identity. T = first time a0 = cos(theta) reaches 0 (the equator, fixed set of the reflection U -> -U).
# Exact: P(T > t) = sum_k c_k exp(-((2k)^2 - 1) t / 2), so the hazard rate tends to 3/2, which is the spectral gap. Compare with the correlator E a0(t) = exp(-3t/2).
import numpy as np
dt = 2e-3; TMAX = 4.0; steps = int(TMAX/dt); N = 20000
TIMES = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
def simulate(bridge, seed=11):
    """Brownian motion on S^3; returns hitting times T, mean a0 (corr) and per-path a0 at the times TIMES (for covariances)."""
    rng = np.random.default_rng(seed)
    q = np.zeros((N, 4)); q[:, 0] = 1; alive = np.ones(N, bool); T = np.full(N, np.inf); corr = np.zeros(steps+1); corr[0] = 1
    snap = {0.0: q[:, 0].copy()}; want = {int(round(t/dt)): t for t in TIMES}
    for s in range(1, steps+1):
        x = q[:, 0].copy(); xi = rng.standard_normal((N, 4))*np.sqrt(dt); xi -= np.sum(xi*q, 1, keepdims=True)*q; q += xi; q /= np.linalg.norm(q, axis=1, keepdims=True)
        y = q[:, 0]; u = rng.uniform(size=N); br = (u < np.exp(-2*np.clip(x*y, 0, None)/dt)) if bridge else False   # uniforms always drawn: common random numbers with and without the bridge; Brownian bridge correction for crossings between steps
        hit = alive & ((y <= 0) | br); T[hit] = s*dt; alive &= ~hit; corr[s] = q[:, 0].mean()
        if s in want: snap[want[s]] = q[:, 0].copy()
    return T, corr, snap, alive
def hazard(T, t0):
    late = T > t0; d = np.isfinite(T[late]).sum(); risk = (np.minimum(T[late], TMAX)-t0).sum(); h = d/risk; return h, h/np.sqrt(d), d
T, corr, snap, alive = simulate(True)
print("paths", N, " still alive at t = 4:", alive.sum())
for t0 in (0.5, 1.0, 1.5):
    h, se, d = hazard(T, t0); print("hazard rate after t0 = %.1f: %.4f +- %.4f (exits %d)   exact gap 1.5" % (t0, h, se, d))
# effective mass between t and t+0.5; delta method with the empirical covariance of a0(t), a0(t+0.5) over the same paths
for t in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5):
    i, j = int(round(t/dt)), int(round((t+0.5)/dt)); ai, aj = snap[t], snap[t+0.5]; ci, cj = ai.mean(), aj.mean(); m = np.log(ci/cj)/0.5 if ci > 0 and cj > 0 else np.nan
    cov = np.cov(ai, aj)  # 2x2
    var = (cov[0, 0]/ci**2 + cov[1, 1]/cj**2 - 2*cov[0, 1]/(ci*cj))/N; se = np.sqrt(var)/0.5
    se_old = np.sqrt((0.25/N)*(1/ci**2+1/cj**2))/0.5
    print("effective mass from E a0 between t = %.1f and %.1f: %.4f +- %.4f   (old independence formula %.4f)" % (t, t+0.5, m, se, se_old))
# bridge bias: paired runs (common random numbers) with and without the correction, over 5 seeds (seed 11 is the run above)
for t0 in (0.5, 1.0, 1.5):
    bias = []
    for sd in (11, 12, 13, 14, 15):
        Tb = simulate(True, sd)[0]; Tn = simulate(False, sd)[0]; hb, hn = hazard(Tb, t0)[0], hazard(Tn, t0)[0]; bias.append(100*(hb-hn)/hb)
        if sd == 11: print("NO bridge, seed 11: hazard rate after t0 = %.1f: %.4f (with bridge %.4f)" % (t0, hn, hb))
    print("bridge bias, hazard after t0 = %.1f: %.1f%% +- %.1f%% (mean and s.e. of 5 seeds; %s)" % (t0, np.mean(bias), np.std(bias, ddof=1)/np.sqrt(5), np.round(bias, 1)))
tt = np.arange(steps+1)*dt; surv = np.array([(T > x).mean() for x in tt[::25]]); np.savez("code/hitting_gap.npz", t=tt[::25], surv=surv, corr=corr[::25], N=N)
