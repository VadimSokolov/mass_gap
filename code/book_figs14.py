# Figure for the chapter on quantum Bayesian computation: the one plaquette gap and the number of levels needed.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from quantum_plaquette import gap, levels
plt.rcParams.update({"font.family": "serif", "font.size": 9, "mathtext.fontset": "cm"})
lam = np.logspace(-1.7, 1, 60); fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7))
ax[0].loglog(lam, [gap(l, 400) for l in lam], "k-", lw=1.2, label="gap"); ax[0].loglog(lam, 0.75*lam, "k--", lw=0.7, label=r"$\frac{3}{4}\lambda$"); ax[0].axhline(np.sqrt(2), color="0.5", ls=":", lw=0.9, label=r"$\sqrt{2}$")
ax[0].set_ylim(0.3, 9); ax[0].set_xlabel(r"$\lambda$"); ax[0].set_ylabel("gap of one plaquette"); ax[0].legend(frameon=False, fontsize=8)
l2 = np.logspace(-1.7, 0.9, 14); ax[1].loglog(l2, [levels(l) for l in l2], "ko", ms=3, label=r"levels for $10^{-8}$"); ax[1].loglog(l2, 6.5/np.sqrt(l2), "k--", lw=0.7, label=r"$6.5\,\lambda^{-1/2}$")
ax[1].set_xlabel(r"$\lambda$"); ax[1].set_ylabel("characters kept"); ax[1].legend(frameon=False, fontsize=8)
plt.tight_layout(); plt.savefig("../figures/quantum_plaquette.pdf")
