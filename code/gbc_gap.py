# GBC for the gap as a parameter. Correlator C(tau) = A exp(-m tau) + B exp(-M tau) + noise, tau = 1..10; prior on (m, A, B, M - m, noise).
# Compare the GBC posterior median of m with the effective mass log C(tau)/C(tau+1) at tau = 3 and tau = 5.
import numpy as np, sys
from sklearn.ensemble import GradientBoostingRegressor as GBR
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 0   # default seed 0 reproduces the original run
rng = np.random.default_rng(5+SEED); tau = np.arange(1, 11)
def sim(N):
    m = rng.uniform(0.2, 1.2, N); A = rng.uniform(0.5, 1.5, N); B = rng.uniform(0, 1.5, N); M = m+rng.uniform(0.5, 2, N); s = 10**rng.uniform(-3, -2, N)
    C = A[:, None]*np.exp(-m[:, None]*tau)+B[:, None]*np.exp(-M[:, None]*tau)+s[:, None]*rng.standard_normal((N, 10)); return m, C
feat = lambda C: np.column_stack([C, C[:, :-1]-C[:, 1:], np.sign(C)*np.log1p(np.abs(C)*1e3)])
if __name__ == "__main__":
    m, C = sim(6000); mc, Cc = sim(2000); mt, Ct = sim(2000)
    g = {q: GBR(loss="quantile", alpha=q, n_estimators=200, max_depth=4, learning_rate=0.08, subsample=0.6, random_state=2+SEED).fit(feat(C), m) for q in (0.05, 0.5, 0.95)}
    lo, md, hi = (g[q].predict(feat(Ct)) for q in (0.05, 0.5, 0.95)); cal = np.quantile(np.maximum(g[0.05].predict(feat(Cc))-mc, mc-g[0.95].predict(feat(Cc))), 0.9)
    print("GBC  median abs error %.4f  rmse %.4f  coverage of calibrated 90%% interval %.3f  mean width %.3f" % (np.median(np.abs(md-mt)), np.sqrt(np.mean((md-mt)**2)), ((mt > lo-cal) & (mt < hi+cal)).mean(), (hi-lo+2*cal).mean()))
    for t in (3, 5):
        with np.errstate(all="ignore"): e = np.log(Ct[:, t-1]/Ct[:, t])
        ok = np.isfinite(e); print("effective mass at tau %d  defined in %.3f of cases  median abs error %.4f  median bias %+.4f" % (t, ok.mean(), np.median(np.abs(e[ok]-mt[ok])), np.median(e[ok]-mt[ok])))
    if SEED == 0: np.savez("code/gbc_gap.npz", mt=mt, md=md, lo=lo-cal, hi=hi+cal, e3=np.log(np.abs(Ct[:, 2]/Ct[:, 3])))
