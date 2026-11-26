"""Covering verification of R(beta) < beta on [0.01, 200] for the heat kernel projected MK map.
If R is increasing and R(beta_{i+1}) < beta_i on a grid, then R(beta) < beta on [beta_i, beta_{i+1}]."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mk_verify import make_R
R = make_R()
g = np.unique(np.concatenate([np.geomspace(0.3, 2, 400), np.linspace(2, 300, 1500)]))
r = np.array([R(b) for b in g])
inc = np.all(np.diff(r) > 0)
cov = np.all(r[1:] < g[:-1])
d = g - r
print("grid points:", len(g), " R increasing on grid:", inc, " covering R(b_{i+1}) < b_i:", cov)
print("min margin b_i - R(b_{i+1}):", (g[:-1]-r[1:]).min())
print("decrement range:", d.min(), d.max(), " argmax beta:", g[np.argmax(d)])
big = g > 20; fit = np.polyfit(1/g[big], d[big], 1); print("fit delta = %.5f + %.4f/beta" % (fit[1], fit[0]))
np.save("/tmp/mk_grid.npy", np.vstack([g, r]))
fig, ax = plt.subplots(figsize=(5.6, 3.5))
ax.semilogx(g, d, "C0", label="heat kernel MK decrement $\\beta-R(\\beta)$")
ax.axhline(0.25, color="k", ls=":", label="$1/4$")
ax.axhline(np.log(2)*11/(3*np.pi**2), color="C3", ls="--", label="one loop $0.2575$")
ax.set_xlabel("$\\beta$ (heat kernel coordinate)"); ax.set_ylim(0, 0.32); ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig("../figures/mk_cover.pdf"); print("saved")
