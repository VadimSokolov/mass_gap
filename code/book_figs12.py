# Figure for the first rung theorem: exact change of the mean count against r for random laws of the class L.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from count_lyapunov import rung
from count_observable import count_law
plt.rcParams.update({"font.family": "serif", "font.size": 9, "mathtext.fontset": "cm"})
rng = np.random.default_rng(5); fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7))
for target, mk, col in ((5, "o", "0.7"), (10, "s", "0.45"), (20, "^", "0.0")):
    R, D = [], []
    for trial in range(60):
        g0 = rng.choice([1.0, 1.5, 3, 8, 25]); gam = g0+rng.exponential(rng.choice([0.01, 1, 4]), 400); n = int(np.searchsorted(np.cumsum(1/(gam+1)), target))+1; gam = gam[:n]
        S1 = (1/(gam+1)).sum(); S2 = (1/(gam+1)**2).sum(); m0, mb, m1, u2 = rung(count_law(gam)); R.append(1-2*S2/S1); D.append(m1-m0)
    ax[0].plot(R, D, mk, ms=3, color=col, mfc="none" if target < 20 else col, label=r"$S_1\approx%d$" % target)
r = np.linspace(0, 1, 50); ax[0].plot(r, -13/16+15/16*r, "k-", lw=0.8); ax[0].axhline(0, color="0.6", lw=0.5); ax[0].set_xlabel(r"$r=1-2S_2/S_1$"); ax[0].set_ylabel("change of mean count, one rung"); ax[0].legend(frameon=False, fontsize=8)
th = np.linspace(0.05, np.pi, 400); ax[1].plot(th, 1/(1-np.cos(th)), "k-", lw=1, label=r"$1/(1-\cos\theta)$"); ax[1].plot(th, np.sin(th)**2/(1-np.cos(th)), "k--", lw=0.8, label="times Haar density")
ax[1].set_xlabel(r"$\theta$"); ax[1].set_ylim(0, 12); ax[1].legend(frameon=False, fontsize=8); ax[1].set_title("the observable whose mean is the count plus two", fontsize=8)
plt.tight_layout(); plt.savefig("../figures/firstrung.pdf")
