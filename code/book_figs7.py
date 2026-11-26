# Colour figures: roots of the heat kernel, ladder coefficients in b, absolute string tension, particle ladder.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
F = "figures/"
G = np.load("/tmp/gammas.npy"); t = G[:, 0]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
for m, c in zip(range(5), cm.plasma(np.linspace(0.05, 0.85, 5))): ax[0].semilogy(t, G[:, m+1], color=c, label=fr"$\gamma_{m+1}$")
for nn, c in zip(range(1, 6), cm.Greys(np.linspace(0.35, 0.8, 5))): ax[0].semilogy(t, np.cosh((2*nn-1)*t/4), color=c, ls=":", lw=1)
ax[0].axhline(1, color="k", lw=0.6); ax[0].set_ylim(0.9, 1e4); ax[0].set_xlabel("$t$"); ax[0].set_ylabel("roots in $-\\cos\\theta$"); ax[0].legend(frameon=False, fontsize=7)
ax[0].set_title(r"roots $\gamma_m(t)$ move outward", fontsize=9)
th = np.linspace(0.01, np.pi-0.01, 300)
for tt, c in zip([0.5, 1, 2, 4], cm.viridis(np.linspace(0.1, 0.9, 4))):
    i = np.argmin(abs(t-tt)); s = sum(1/(g+np.cos(th)) for g in G[i, 1:]); ax[1].plot(th, -np.sin(th)*s, color=c, label=f"$t={t[i]:.1f}$")
ax[1].set_xlabel(r"$\theta$"); ax[1].set_ylabel(r"$\partial_\theta\log K_t$, first five roots"); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title("score increases with $t$", fontsize=9)
fig.tight_layout(); fig.savefig(F+"heat_roots.pdf"); plt.close()
L = np.load("/tmp/ladder_b.npy"); fig, ax = plt.subplots(figsize=(5.2, 3.4))
sc = ax.scatter(L[:, 1], L[:, 2], c=L[:, 0], cmap="turbo", s=45, zorder=3); ax.plot(L[:, 1], L[:, 2], color="#7f8c8d", lw=1)
for b, a, c, _ in L:
    if b in (1.2, 1.6, 2.0, 3.0): ax.annotate(f"$b={b:g}$", (a, c), textcoords="offset points", xytext=(6, 4), fontsize=7)
ax.plot([0.2575], [0.0403], "*", color="#c0392b", ms=14, label="perturbation theory"); ax.legend(frameon=False, fontsize=8)
ax.set_xlabel("one loop coefficient $a(b)$"); ax.set_ylabel("two loop coefficient $c(b)$"); fig.colorbar(sc, label="scale factor $b$")
fig.tight_layout(); fig.savefig(F+"ladder_b.pdf"); plt.close()
s4 = np.load("/tmp/mk_sigma4.npy"); s3 = np.load("/tmp/mk_sigma3.npy")
b4 = np.array([2.2, 2.3, 2.4, 2.5, 2.635]); d4 = np.array([0.4690, 0.3690, 0.2660, 0.1871, 0.1207])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
ax[0].semilogy(s4[:, 0], s4[:, 1], color="#1f4e79", label="ladder"); ax[0].semilogy(s4[:, 0]+1.45, s4[:, 1], color="#1f4e79", ls="--", label=r"ladder, $\beta$ shifted by 1.45")
ax[0].semilogy(b4, d4, "o", color="#c0392b", label="lattice"); xs = np.linspace(0.5, 2.2, 50); ax[0].semilogy(xs, np.sqrt(-np.log(xs/4)), color="#16a085", ls=":", label="strong coupling")
ax[0].set_xlabel(r"$\beta$"); ax[0].set_ylabel(r"$a\sqrt{\sigma}$"); ax[0].set_title("four dimensions", fontsize=9); ax[0].legend(frameon=False, fontsize=6.5); ax[0].set_xlim(0.5, 4.2)
ax[1].plot(s3[:, 0], s3[:, 0]*s3[:, 1], "o-", ms=3, color="#8e44ad", label="ladder"); ax[1].axhline(1.341, color="#c0392b", label="lattice continuum")
ax[1].set_ylim(0, 1.6); ax[1].set_xlabel(r"$\beta$"); ax[1].set_ylabel(r"$\beta\,a\sqrt{\sigma}$"); ax[1].set_title("three dimensions", fontsize=9); ax[1].legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(F+"mk_sigma.pdf"); plt.close()
H = np.load("/tmp/pmk_hist.npy"); T = np.load("/tmp/pmk_tr.npy"); fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
for s, c in zip(range(6), cm.viridis(np.linspace(0.05, 0.92, 6))): ax[0].hist(H[s], bins=80, density=True, histtype="step", color=c, label=f"step {s}")
ax[0].plot(th, (2/np.pi)*np.sin(th)**2, "k--", lw=1); ax[0].set_xlabel(r"$\theta$"); ax[0].set_ylabel("density of the cloud"); ax[0].legend(frameon=False, fontsize=6.5)
for r, c, l in zip(T, ["#1f4e79", "#c0392b"], [r"$\beta=3$", r"$\beta=2$"]):
    ax[1].plot(r[:6], "o", color=c, label="particles, "+l); ax[1].plot(r[6:], "-", color=c, label="exact, "+l)
ax[1].set_xlabel("ladder step"); ax[1].set_ylabel("$u_2$"); ax[1].legend(frameon=False, fontsize=6.5)
fig.tight_layout(); fig.savefig(F+"particle_mk.pdf"); plt.close(); print("ok")
