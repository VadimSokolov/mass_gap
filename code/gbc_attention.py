# Attention inside GBC. Data: a set of n = 256 plaquettes a_i (two dimensional SU(2)); target: the posterior of beta.
# (a) Pooling: tokens are embedded by the characters chi_2..chi_9 of the plaquette and mean pooled (attention with uniform weights);
#     a quantile head is trained on the pooled vector. The sufficient statistic is the first pooled coordinate; the network must find it.
# (b) In context attention: the query's pooled vector attends to the pooled vectors of the training sets with weights
#     softmax(-|W(s - s_i)|^2 / 2), W diagonal and learned on a validation split; posterior quantiles are attention weighted quantiles
#     of the training couplings. This is a one layer, one head transformer acting as a kernel smoother.
import numpy as np, sys
from sklearn.ensemble import GradientBoostingRegressor as GBR
sys.path.insert(0, "code"); import gbc_exact as G
qs = G.qs; n = G.n
def embed(a): th = np.arccos(np.clip(a, -1, 1)); th = np.clip(th, 1e-9, np.pi-1e-9); return np.array([(np.sin(k*th)/np.sin(th)).mean()/k for k in range(2, 10)])
def sim(N): b = G.rng.uniform(0.5, 6, N); return b, np.array([embed(G.draw(x, n)) for x in b])
def wq(w, v, q): o = np.argsort(v); c = np.cumsum(w[o]); return np.interp(q, c/c[-1], v[o])
def head(Xtr, b, Xte): return np.column_stack([GBR(loss="quantile", alpha=q, n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8, random_state=1).fit(Xtr, b).predict(Xte) for q in qs])
def run(seed, verbose):
    """One configuration: 5000 training sets, 300 validation, 500 test sets, all generators on the same sets."""
    G.rng = np.random.default_rng(seed)
    b, S = sim(5000); bv, Sv = sim(300); bt, St = sim(500); ex = np.array([G.exact_q(s[0]) for s in St]); sd = S.std(0)
    P = head(S, b, St); Ps = head(S[:, [0]], b, St[:, [0]])
    def attend(Sq, scale):
        out = []
        for s in Sq: d = ((s-S)/(sd*scale))**2; w = np.exp(-0.5*d.sum(1)); out.append([wq(w, b, q) for q in qs])
        return np.array(out)
    exv = np.array([G.exact_q(s[0]) for s in Sv]); best = None
    for first in (0.03, 0.06, 0.12, 0.25):
        for rest in (0.5, 2.0, 1e6):
            sc = np.array([first]+[rest]*7); e = np.abs(attend(Sv, sc)-exv).mean()
            if best is None or e < best[0]: best = (e, first, rest)
    A = attend(St, np.array([best[1]]+[best[2]]*7))
    rows = {}
    for name, Q in (("sufficient statistic head", Ps), ("pooled characters, boosted head", P), ("pooled characters, attention head", A), ("exact posterior", ex)):
        err = np.abs(Q-ex).mean(0); rows[name] = (err, (Q[:, 4]-Q[:, 0]).mean(), ((bt > Q[:, 0]) & (bt < Q[:, 4])).mean())
        if verbose: print("seed %d  %-34s quantile error %s  90%% width %.3f  coverage %.3f" % ((seed, name, np.round(err, 3)) + rows[name][1:]))
    if verbose: print("seed %d  learned attention scales: first coordinate %.2f, others %g" % (seed, best[1], best[2]))
    return rows
if __name__ == "__main__":
    SEEDS = (12, 13, 14, 15, 16); R = [run(sd, True) for sd in SEEDS]   # seed 12 is the configuration used before the replicates
    print("\nmean over %d replicate seeds (s.e. in brackets): quantile error averaged over the five quantiles, 90%% width, coverage" % len(SEEDS))
    for name in R[0]:
        e = np.array([r[name][0].mean() for r in R]); w = np.array([r[name][1] for r in R]); c = np.array([r[name][2] for r in R]); f = lambda x: "%.3f (%.3f)" % (x.mean(), x.std(ddof=1)/np.sqrt(len(x)))
        er = np.array([r[name][0] for r in R]).mean(0)
        print("%-34s error %s  [per quantile, mean over seeds %s]  width %s  coverage %s" % (name, f(e), np.round(er, 3), f(w), f(c)))
