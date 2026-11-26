"""Vertical likelihood summary, implied likelihood ordinate distributions, MK phase portrait, roadmap. Polson and Sokolov."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
out = "../figures/"
# errors by method and seed: all three seeds (5, 11, 23) from vl_methods.py, saved by it in ../results/vl_res_<seed>.npy
meth = ["prior IS", "harmonic mean", "power posterior", "nested", "vertical likelihood", "split sampling"]
keys = ["1 importance sampling (prior)", "4 harmonic mean (w = 1)", "3 power posterior / path sampling", "7 nested sampling",
        "8 vertical likelihood (weighted slice)", "9 split sampling (Birge Chang Polson)"]
E = []
for sd in (5, 11, 23):
    r = np.load("../results/vl_res_%d.npy" % sd, allow_pickle=True).item(); E.append([r[k][1] for k in keys])
E = np.array(E); rmse = np.sqrt(np.mean(E**2, 0)); print(dict(zip(meth, np.round(rmse, 3))))
np.save("../results/vl_rmse.npy", np.vstack([E, rmse]))
fig, ax = plt.subplots(figsize=(6.2, 3.2)); x = np.arange(len(meth)); cols = ["#4C72B0", "#DD8452", "#55A868"]
for k in range(3): ax.bar(x+(k-1)*0.25, E[k], width=0.25, color=cols[k], label="seed %d" % [5, 11, 23][k])
ax.plot(x, rmse, "k_", ms=18, mew=2, label="RMSE")
ax.axhline(0, color="k", lw=0.6); ax.set_xticks(x); ax.set_xticklabels(meth, rotation=25, ha="right", fontsize=8)
ax.set_ylabel("error in $\\log Z(2)$"); ax.legend(frameon=False, fontsize=7, ncol=4)
fig.tight_layout(); fig.savefig(out+"vl_errors.pdf"); plt.close()
# implied distributions of the action under prior, posterior and the default vertical likelihood weight, from a nested run
import vl_core as v
rng = v.rng; L, P = v.L, v.P; N = 100; live = [v.haar((L, L, 2)) for _ in range(N)]; Sl = np.array([v.action_all(u) for u in live]); dS, lX = [], []; eps = 0.6
for it in range(3000):
    i = np.argmax(Sl); Smax = Sl[i]; dS.append(Smax); lX.append(-(it+1)/N)
    j = rng.integers(N-1); j = j if j < i else j+1; U = live[j].copy(); Sj = Sl[j]; acc = 0
    for _ in range(50):
        x_, y_, mu = rng.integers(L), rng.integers(L), rng.integers(2); old = v.local_S(U, x_, y_, mu); keep = U[x_, y_, mu].copy()
        r = rng.standard_normal(3)*eps; th = np.linalg.norm(r); U[x_, y_, mu] = v.qmul(np.concatenate([[np.cos(th)], np.sin(th)*r/max(th, 1e-12)]), keep)
        Sn = Sj-old+v.local_S(U, x_, y_, mu)
        if Sn < Smax: Sj = Sn; acc += 1
        else: U[x_, y_, mu] = keep
    eps = min(max(eps*np.exp(0.5*(acc/50-0.35)), 0.02), 1.5); live[i] = U; Sl[i] = Sj
dS = np.array(dS); lX = np.array(lX); X = np.exp(lX); dX = -np.diff(np.concatenate([[1.0], X]))
fig, ax = plt.subplots(1, 2, figsize=(6.4, 3.0))
ax[0].plot(dS, lX, "C0"); ax[0].set_xlabel("action level $s$"); ax[0].set_ylabel("$\\log X(s)$, prior mass of $\\{S<s\\}$")
eta = np.exp(-30); Wd = 1/np.maximum(eta, X)
for wts, lab, c in [(dX, "prior", "C0"), (dX*np.exp(-2*dS), "posterior, $\\beta=2$", "C3"), (dX*Wd, "vertical likelihood $W=1/\\max(\\eta,X)$", "C2")]:
    h, edges = np.histogram(dS, bins=40, weights=wts, density=True); ax[1].step(0.5*(edges[1:]+edges[:-1]), h, color=c, label=lab)
ax[1].set_xlabel("action $S$"); ax[1].set_ylabel("implied density"); ax[1].legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(out+"vl_implied.pdf"); plt.close()
# MK phase portrait in (u2, u3)
from mk_flow import ive as _i
import mk_flow as mk
fig, ax = plt.subplots(figsize=(5.2, 3.6)); cm = plt.cm.viridis
for k, b0 in enumerate([0.8, 1.5, 3, 5, 8, 12]):
    u = mk.ive(mk.n, b0)/mk.ive(1, b0); pts = [(u[1], u[2])]
    for s in range(60):
        u = mk.step(u); pts.append((u[1], u[2]))
        if u[1] < 1e-3: break
    pts = np.array(pts); ax.plot(pts[:, 0], pts[:, 1], "-o", ms=2.5, color=cm(k/6), label="$\\beta_0=%g$" % b0)
ub = np.linspace(0.01, 0.99, 200); ax.plot(ub, ub**(8/3), "k:", lw=0.8, label="heat kernel family")
ax.set_xlabel("$u_2$"); ax.set_ylabel("$u_3$"); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(out+"mk_phase.pdf"); plt.close()
# roadmap diagram
fig, ax = plt.subplots(figsize=(6.4, 3.6)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
boxes = [(0.2, 4.3, "(i) positivity\nat every scale", "#55A868", "done"), (2.2, 4.3, "(ii) strong\ncoupling", "#55A868", "done"),
         (4.2, 4.3, "(iii) MK ladder\nand decrement", "#C9B33F", "nearly done"), (6.2, 4.3, "(iv) MK versus $\\mathbb{Z}^4$\nsandwich", "#C44E52", "open"),
         (8.2, 4.3, "(v) continuum\nand axioms", "#4C72B0", "available"),
         (1.2, 1.4, "dynamic gap\n(log Sobolev,\nrapid mixing)", "#8172B2", "route"), (4.2, 1.4, "latent variable\ncoupling\n(Edwards Sokal)", "#8172B2", "route"),
         (7.2, 1.4, "Bayesian gap\nestimation\n(lattice data)", "#937860", "diagnostic")]
for x, y, t, c, st in boxes:
    ax.add_patch(FancyBboxPatch((x, y), 1.6, 1.2, boxstyle="round,pad=0.05", fc=c, ec="none", alpha=0.85))
    ax.text(x+0.8, y+0.68, t, ha="center", va="center", fontsize=7.5, color="white"); ax.text(x+0.8, y+0.15, st, ha="center", fontsize=7, color="white", style="italic")
for x in [1.8, 3.8, 5.8, 7.8]: ax.annotate("", xy=(x+0.4, 4.9), xytext=(x, 4.9), arrowprops=dict(arrowstyle="->", color="k"))
for x0 in [2.0, 5.0]: ax.annotate("", xy=(7.0, 4.25), xytext=(x0, 2.65), arrowprops=dict(arrowstyle="->", color="#8172B2", lw=1.2))
ax.annotate("", xy=(7.0, 4.25), xytext=(8.0, 2.65), arrowprops=dict(arrowstyle="->", color="#937860", lw=1.0, ls="--"))
fig.tight_layout(); fig.savefig(out+"roadmap.pdf"); plt.close()
print("ok")
