"""All vertical likelihood methods (Polson and Scott 2015) plus split sampling (Birge, Chang and Polson)
on 2D SU(2) lattice gauge theory, 3x3 torus; target log Z(2) = log E_Haar exp(-2 S), exact by characters.
Each method gets roughly the same budget of single link updates. Polson and Sokolov, Nov 2026."""
import numpy as np, time
from scipy.special import logsumexp
import vl_core as v, sys
v.rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
rng = v.rng; L, P = v.L, v.P; B0 = 2.0; EX = v.logZ_exact(B0); BUDGET = 240000
def sweep_metro(U, S, logw, nprop, eps=0.9):
    """Metropolis on single links with target Haar x exp(logw(S))"""
    for _ in range(nprop):
        x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
        old = v.local_S(U, x, y, mu); keep = U[x, y, mu].copy()
        r = rng.standard_normal(3)*eps; th = np.linalg.norm(r)
        U[x, y, mu] = v.qmul(np.concatenate([[np.cos(th)], np.sin(th)*r/max(th, 1e-12)]), keep)
        Sn = S - old + v.local_S(U, x, y, mu)
        if np.log(rng.random()) < logw(Sn) - logw(S): S = Sn
        else: U[x, y, mu] = keep
    return U, S
res = {}
def report(name, est, t0, extra=""):
    res[name] = (est, est-EX, time.time()-t0); print("%-34s logZ=%9.4f  error=%+.4f  %4.0fs %s" % (name, est, est-EX, time.time()-t0, extra))
print("exact log Z(2) = %.4f" % EX)
# 1 importance sampling from the prior
t0 = time.time(); Sp = np.array([v.action_all(v.haar((L, L, 2))) for _ in range(BUDGET//18)])
report("1 importance sampling (prior)", logsumexp(-B0*Sp)-np.log(len(Sp)), t0)
# 2 slice sampling heat bath: posterior draws, used by 4
t0 = time.time(); U = v.slice_heatbath(v.haar((L, L, 2)), B0, 50); post = []
for i in range(BUDGET//18//2): U = v.slice_heatbath(U, B0, 1); post.append(v.action_all(U))
post = np.array(post); print("2 slice heat bath: posterior mean S %.3f (exact %.3f), %.0fs" % (post.mean(), v.meanS_exact(B0), time.time()-t0))
# 4 harmonic mean (w = 1)
report("4 harmonic mean (w = 1)", -(logsumexp(B0*post)-np.log(len(post))), time.time())
# 3 annealing / power posterior thermodynamic integration (w(u) = a u^(a-1))
t0 = time.time(); rungs = B0*np.linspace(0, 1, 11)**2; ms = []
for b in rungs:
    if b == 0: ms.append(np.mean(Sp)); continue
    U = v.slice_heatbath(v.haar((L, L, 2)), b, 20); ss = []
    for i in range(BUDGET//18//2//10): U = v.slice_heatbath(U, b, 1); ss.append(v.action_all(U))
    ms.append(np.mean(ss))
report("3 power posterior / path sampling", -np.trapezoid(ms, rungs), t0)
# 6 energy level: Wang Landau (multicanonical, w = 1/N(s)) then production
t0 = time.time(); edges = np.linspace(0.5, 15.5, 41); nb = len(edges)-1; lg = np.zeros(nb); H = np.zeros(nb); f = 1.0
bin_ = lambda s: min(max(int(np.searchsorted(edges, s)-1), 0), nb-1)
U = v.haar((L, L, 2)); S = v.action_all(U)
while not (edges[0] <= S < edges[-1]): U = v.haar((L, L, 2)); S = v.action_all(U)
used = 0
while f > 2e-3 and used < BUDGET:
    for _ in range(3000):
        x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
        old = v.local_S(U, x, y, mu); keep = U[x, y, mu].copy(); U[x, y, mu] = v.haar(())
        Sn = S - old + v.local_S(U, x, y, mu); i0, i1 = bin_(S), bin_(Sn)
        if edges[0] <= Sn < edges[-1] and np.log(rng.random()) < lg[i0]-lg[i1]: S = Sn; i0 = i1
        else: U[x, y, mu] = keep
        lg[i0] += f; H[i0] += 1
    used += 3000; vis = H > 0
    if H.min() > 0.6*H.mean(): f /= 2; H[:] = 0
vis = lg > 0; Sc = 0.5*(edges[1:]+edges[:-1])
# anchor the density of states to the prior: mass of the top bins from the prior sample
m_hi = np.mean(Sp > edges[-2]); lgn = lg - lg[-1] + np.log(max(m_hi, 1e-12))
report("6 energy levels (Wang Landau)", logsumexp(lgn[vis]-B0*Sc[vis]), t0, "(%d updates, final f=%.3g)" % (used, f))
# 7 nested sampling
t0 = time.time(); N = 100; live = [v.haar((L, L, 2)) for _ in range(N)]; Sl = np.array([v.action_all(u) for u in live])
dS, lX = [], []; it = 0; eps = 0.6; nprop = 50
while True:
    i = np.argmax(Sl); Smax = Sl[i]; dS.append(Smax); lX.append(-(it+1)/N)
    j = rng.integers(N-1); j = j if j < i else j+1; U = live[j].copy(); Sj = Sl[j]; acc = 0
    for _ in range(nprop):
        x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
        old = v.local_S(U, x, y, mu); keep = U[x, y, mu].copy(); r = rng.standard_normal(3)*eps; th = np.linalg.norm(r)
        U[x, y, mu] = v.qmul(np.concatenate([[np.cos(th)], np.sin(th)*r/max(th, 1e-12)]), keep)
        Sn = Sj - old + v.local_S(U, x, y, mu)
        if Sn < Smax: Sj = Sn; acc += 1
        else: U[x, y, mu] = keep
    eps = min(max(eps*np.exp(0.5*(acc/nprop-0.35)), 0.02), 1.5); live[i] = U; Sl[i] = Sj; it += 1
    if it*nprop > BUDGET: break
dS = np.array(dS); lX = np.array(lX)
lw = np.log(np.exp(np.concatenate([[0], lX[:-1]]))-np.exp(lX))
allS = np.concatenate([dS, np.sort(Sl)[::-1]]); allw = np.concatenate([lw, np.full(N, lX[-1]-np.log(N))])
report("7 nested sampling", logsumexp(allw-B0*allS), t0, "(N=%d, %d iterations)" % (N, it))
# trapezoid in X (Yakowitz style) versus rectangle for the same nested run
Xs = np.exp(np.concatenate([[0], lX])); Ls = np.exp(-B0*dS)
trap = np.sum(0.5*(np.concatenate([[np.exp(-B0*2*P)], Ls[:-1]]) + Ls)*(Xs[:-1]-Xs[1:])) + np.exp(lX[-1])*np.mean(np.exp(-B0*Sl))
print("   nested quadrature: rectangle %.4f  trapezoid %.4f" % (logsumexp(allw-B0*allS), np.log(trap)))
Xhat = lambda s: np.exp(np.interp(-s, -dS, lX)) if s < dS[0] else 1.0          # pilot prior mass of {S < s}
# 8 vertical likelihood, weighted slice sampling with the default weight W = 1/max(eta, Z(u))
t0 = time.time(); eta = np.exp(lX[-1]); logW = lambda s: -np.log(max(eta, Xhat(s)))
U = v.haar((L, L, 2)); S = v.action_all(U); U, S = sweep_metro(U, S, logW, 20000); Ss = []
for i in range((BUDGET-20000)//18): U, S = sweep_metro(U, S, logW, 18); Ss.append(S)
Ss = np.array(Ss); lw8 = np.array([logW(s) for s in Ss])
report("8 vertical likelihood (weighted slice)", logsumexp(-B0*Ss-lw8)-logsumexp(-lw8), t0, "(eta = e^%.1f, pilot from 7)" % np.log(eta))
# 9 split sampling (Birge, Chang and Polson): levels s_t, joint (x, t), adaptive weights omega_t
t0 = time.time(); T = 13; lev = np.concatenate([[2*P+1e-9], np.interp(np.linspace(-1.2, -14, T-1), lX[::-1], dS[::-1])]); T = len(lev); lom = np.zeros(T); Hs = np.zeros(T); fs = 1.0
U = v.haar((L, L, 2)); S = v.action_all(U); t = 0; band = [[] for _ in range(T)]; used = 0
while used < BUDGET:
    for _ in range(2000):
        x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
        old = v.local_S(U, x, y, mu); keep = U[x, y, mu].copy(); r = rng.standard_normal(3)*min(1.2, 0.15+0.9*np.sqrt(lev[t]/(2*P))); th = np.linalg.norm(r)
        U[x, y, mu] = v.qmul(np.concatenate([[np.cos(th)], np.sin(th)*r/max(th, 1e-12)]), keep)
        Sn = S - old + v.local_S(U, x, y, mu)
        if Sn < lev[t]: S = Sn
        else: U[x, y, mu] = keep
        tp = t + (1 if rng.random() < 0.5 else -1)
        if 0 <= tp < T and S < lev[tp] and np.log(rng.random()) < lom[tp]-lom[t]: t = tp
        lom[t] -= fs; Hs[t] += 1
        if fs < 0.13: band[t].append(S)
    used += 2000
    if Hs.min() > 0.6*Hs.mean(): fs /= 2; Hs[:] = 0
# occupancy flat => omega_t proportional to 1/X(s_t); X(s_0) = 1
lX9 = -(lom - lom[0])
Xl = np.exp(lX9); Z9 = 0.0
for k in range(T):
    Xnext = Xl[k+1] if k+1 < T else 0.0
    sk = np.array(band[k]); sk = sk[sk >= (lev[k+1] if k+1 < T else -1)] if len(sk) else sk
    m = np.mean(np.exp(-B0*sk)) if len(sk) else np.exp(-B0*lev[k])
    Z9 += (Xl[k]-Xnext)*m
report("9 split sampling (Birge Chang Polson)", np.log(Z9), t0, "(%d levels, final f=%.3g)" % (T, fs))
import os; os.makedirs("../results", exist_ok=True)
np.save("../results/vl_res_%s.npy" % (sys.argv[1] if len(sys.argv) > 1 else "5"), res)
