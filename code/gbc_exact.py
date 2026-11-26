# GBC against the exact posterior. Data: n = 256 independent SU(2) plaquettes with density sqrt(1-a^2) exp(beta a); prior beta ~ U(0.5, 6).
# Quantile generators (gradient boosting) trained on the sufficient statistic (mean of a) and on a lossy one bit statistic (fraction a > 1/2).
import numpy as np
from scipy.special import ive
from sklearn.ensemble import GradientBoostingRegressor as GBR
rng = np.random.default_rng(3); n = 256; grid = np.linspace(0.5, 6, 1101); qs = (0.05, 0.25, 0.5, 0.75, 0.95)
def draw(beta, n):
    out = np.empty(0)
    while out.size < n:
        a = rng.uniform(-1, 1, 4*n); keep = rng.uniform(size=4*n) < np.sqrt(1-a*a)*np.exp(beta*(a-1)); out = np.concatenate([out, a[keep]])
    return out[:n]
def exact_q(abar):
    ll = n*(grid*abar-np.log(ive(1, grid)/grid)-grid); w = np.exp(ll-ll.max()); c = np.cumsum(w)/w.sum(); return np.interp(qs, c, grid)
def sim(N):
    b = rng.uniform(0.5, 6, N); S = np.array([(lambda a: (a.mean(), (a > 0.5).mean()))(draw(x, n)) for x in b]); return b, S
if __name__ == "__main__":
    b, S = sim(6000); bt, St = sim(600); ex = np.array([exact_q(s) for s in St[:, 0]]); res = {}
    print("exact posterior: coverage of the exact 90%% interval %.3f (se %.3f, %d test sets), width %.3f" % (((bt > ex[:, 0]) & (bt < ex[:, 4])).mean(), np.sqrt(0.09/len(bt)), len(bt), (ex[:, 4]-ex[:, 0]).mean()))
    for name, col in (("sufficient", 0), ("lossy", 1)):
        P = np.column_stack([GBR(loss="quantile", alpha=q, n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8, random_state=1).fit(S[:, [col]], b).predict(St[:, [col]]) for q in qs])
        res[name] = P
        print("%-10s mean abs quantile error vs exact %s  width of 90%% interval %.3f (exact %.3f)  coverage %.3f" % (name, np.round(np.abs(P-ex).mean(0), 3), (P[:, 4]-P[:, 0]).mean(), (ex[:, 4]-ex[:, 0]).mean(), ((bt > P[:, 0]) & (bt < P[:, 4])).mean()))
    np.savez("code/gbc_exact.npz", bt=bt, St=St, ex=ex, suff=res["sufficient"], lossy=res["lossy"])
