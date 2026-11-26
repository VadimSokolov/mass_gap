# Figure: the negative binomial estimator M_k of the gap and the effective mass, one plaquette at lambda = 1/2.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import window_gap as W
plt.rcParams.update({"font.family": "serif", "font.size": 9, "mathtext.fontset": "cm"})
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7)); ks = np.arange(1, 41)
for sigma, ls in ((1.0, "-"), (3.0, "--"), (8.0, ":")): ax[0].semilogy(ks, [W.nb(W.C, k, sigma)-W.gap for k in ks], "k"+ls, lw=1, label=r"$q=e^{-%ga}$" % sigma)
ax[0].set_xlabel(r"$k$"); ax[0].set_ylabel(r"$M_k(q)-E_0$"); ax[0].legend(frameon=False, fontsize=8)
rng = np.random.default_rng(3); noise = 1e-4*W.C[0]; t = np.arange(0.25, 5.01, 0.25)
for r in range(12): ax[1].plot(t, [W.eff(W.C+noise*rng.standard_normal(W.C.size), s) for s in t], color="0.7", lw=0.5)
ks2 = np.arange(2, 30); cen = [k*W.a/(np.exp(3*W.a+W.gap*W.a)-1) for k in ks2]
for r in range(12): Cn = W.C+noise*rng.standard_normal(W.C.size); ax[1].plot(cen, [W.nb(Cn, k, 3.0) for k in ks2], color="k", lw=0.5)
ax[1].axhline(W.gap, color="k", ls=":", lw=0.8); ax[1].set_ylim(1.21, 1.27); ax[1].set_xlim(0, 5); ax[1].set_xlabel(r"$\tau$, or the centre of the kernel"); ax[1].set_ylabel("estimate of the gap")
plt.tight_layout(); plt.savefig("../figures/window_gap.pdf")
