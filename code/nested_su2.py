"""Nested sampling (Skilling) for 2D SU(2) lattice gauge theory on an L x L torus, Haar prior,
likelihood exp(-beta S), S = sum_p (1 - tr U_p / 2). One run gives Z(beta) for all beta.
Exact: Z(beta) = exp(-beta P) sum_n c_n(beta)^P with c_n = 2 I_n(beta)/beta (torus, chi = 0).
Polson and Sokolov, Nov 2026."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ive, logsumexp
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(11); L = 3; P = L*L
def qmul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], -1)
def conj(a): return a*np.array([1, -1, -1, -1])
def haar(shape):
    v = rng.standard_normal(shape+(4,)); return v/np.linalg.norm(v, axis=-1, keepdims=True)
def plaq(U, x, y):   # U[x,y,mu,4]
    a = U[x, y, 0]; b = U[(x+1) % L, y, 1]; c = conj(U[x, (y+1) % L, 0]); d = conj(U[x, y, 1])
    return qmul(qmul(qmul(a, b), c), d)[..., 0]          # = tr U_p / 2
def action(U):
    return sum(1-plaq(U, x, y) for x in range(L) for y in range(L))
def plaqs_of_link(x, y, mu):
    return [(x, y), (x, (y-1) % L)] if mu == 0 else [(x, y), ((x-1) % L, y)]
def small_rot(eps):
    v = rng.standard_normal(3)*eps; th = np.linalg.norm(v)
    return np.concatenate([[np.cos(th)], np.sin(th)*v/max(th, 1e-12)])
def constrained_walk(U, S, Smax, nprop, eps):
    acc = 0
    for _ in range(nprop):
        x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
        ps = plaqs_of_link(x, y, mu); old = sum(1-plaq(U, *p) for p in ps)
        keep = U[x, y, mu].copy(); U[x, y, mu] = qmul(small_rot(eps), keep)
        new = sum(1-plaq(U, *p) for p in ps); Snew = S - old + new
        if Snew < Smax: S = Snew; acc += 1
        else: U[x, y, mu] = keep
    return U, S, acc/nprop
N = 200; live = [haar((L, L, 2)) for _ in range(N)]; Sl = np.array([action(U) for U in live])
dead_S, logX = [], []; eps = 0.6; it = 0
while True:
    i = np.argmax(Sl); Smax = Sl[i]; dead_S.append(Smax); logX.append(-(it+1)/N)
    j = rng.integers(N-1); j = j if j < i else j+1
    U = live[j].copy(); U, S, a = constrained_walk(U, Sl[j], Smax, 160, eps)
    eps *= np.exp(0.5*(a-0.35)); eps = min(max(eps, 0.02), 1.5)
    live[i] = U; Sl[i] = S; it += 1
    if Smax < 0.6 or it > 12000: break
dead_S = np.array(dead_S); logX = np.array(logX)
# add the remaining live points
logw = np.log(np.exp(np.concatenate([[0], logX[:-1]])) - np.exp(logX))
liveS = np.sort(Sl)[::-1]; logw_live = np.full(N, logX[-1]-np.log(N))
allS = np.concatenate([dead_S, liveS]); allw = np.concatenate([logw, logw_live])
betas = np.linspace(0.05, 8, 80)
logZ_ns = np.array([logsumexp(allw - b*allS) for b in betas])
def logZ_exact(b):
    n = np.arange(1, 400); logc = np.log(2*ive(n, b)/b) + b     # log c_n = log(2 I_n/b)
    return -b*P + logsumexp(P*logc)
logZ_ex = np.array([logZ_exact(b) for b in betas])
print("iterations", it, " final Smax", Smax, " max |logZ_ns - exact| =", np.max(np.abs(logZ_ns-logZ_ex)))
for b in [0.5, 1, 2, 4, 8]:
    k = np.argmin(abs(betas-b)); print("beta=%.2f  NS %.4f  exact %.4f" % (betas[k], logZ_ns[k], logZ_ex[k]))
# plain Monte Carlo cross check at beta = 1
mc = np.mean([np.exp(-action(haar((L, L, 2)))) for _ in range(20000)]); print("plain MC log Z(1) =", np.log(mc), " exact", logZ_exact(1.0))
# mean action <S>_beta from NS vs exact derivative
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.9))
ax[0].plot(betas, logZ_ex, "k", lw=2.5, alpha=0.35, label="exact")
ax[0].plot(betas, logZ_ns, "C0--", label="nested sampling")
ax[0].set_xlabel("$\\beta$"); ax[0].set_ylabel("$\\log Z(\\beta)$"); ax[0].legend(frameon=False, fontsize=8)
ax[1].plot(-logX, dead_S, "C3", lw=0.8)
ax[1].set_xlabel("$-\\log X$ (prior mass compressed)"); ax[1].set_ylabel("action threshold $S^*$")
fig.tight_layout(); fig.savefig("../figures/nested_su2.pdf"); print("ok")
