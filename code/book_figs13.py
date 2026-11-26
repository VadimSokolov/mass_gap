# Figure: the Thorin measure of the Hartman Watson law, and its distribution function against the zeros of K_{iy}.
import numpy as np, mpmath as mp, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hw_thorin import dtheta, phase
plt.rcParams.update({"font.family": "serif", "font.size": 9, "mathtext.fontset": "cm"})
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7)); y = np.linspace(0.02, 12, 240)
for r, ls in ((0.5, "-"), (2.0, "--"), (8.0, ":")): ax[0].plot(y, [-dtheta(v, r)/np.pi for v in y], "k"+ls, lw=1, label=r"$r=%g$" % r)
ax[0].plot(y, [float(mp.euler+mp.re(mp.digamma(mp.mpc(1, v))))/np.pi for v in y], color="0.6", lw=0.8, label="Bondesson's law"); ax[0].set_xlabel(r"$y=\sqrt{2t}$"); ax[0].set_ylabel(r"Thorin density $-\theta'(y)/\pi$"); ax[0].legend(frameon=False, fontsize=8)
ys = np.linspace(0, 14, 2801); th = phase(1.0, ys); ax[1].plot(ys, -th/np.pi, "k-", lw=1)
kz = [2.9625, 4.5345, 5.8799, 7.1076, 8.2589, 9.3551, 10.4087, 11.4281, 12.419, 13.3859]
for n, z in enumerate(kz): ax[1].plot([0, z, z], [n+1, n+1, 0], color="0.7", lw=0.5)
ax[1].plot(kz, np.arange(1, 11), "ko", ms=3, mfc="w"); ax[1].set_xlim(0, 14); ax[1].set_ylim(0, 10.8); ax[1].set_xlabel(r"$y$"); ax[1].set_ylabel(r"$U([0,y^2/2])=-\theta(y)/\pi$, $r=1$")
plt.tight_layout(); plt.savefig("../figures/hw_thorin.pdf")
