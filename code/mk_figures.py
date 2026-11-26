"""MK flow figures. Polson and Sokolov, Nov 2026."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mk_flow import flow
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
cAF = np.log(2)*11/(3*np.pi**2)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.7))
B, D = [], []
for b0 in [1.5, 3, 6, 12, 24, 40]:
    f = flow(b0, steps=400, stop=0.02)
    bb, dd = f[1:-1], f[1:-1]-f[2:]
    B += list(bb); D += list(dd)
    ax[0].plot(bb, dd, ".", ms=3, color="C0")
B, D = np.array(B), np.array(D)
ax[0].axhline(0.25, color="k", ls=":", label="MK weak coupling: $1/4$")
ax[0].axhline(cAF, color="C3", ls="--", label="one loop AF: $11\\ln 2/(3\\pi^2)$")
ax[0].set_xscale("log"); ax[0].set_xlabel("$\\beta_{\\mathrm{eff}}$"); ax[0].set_ylabel("decrement $\\delta$ per step $a\\to2a$")
ax[0].set_ylim(0, 0.6); ax[0].legend(frameon=False, fontsize=9)
b0s = np.linspace(0.5, 30, 45); N = np.array([len(flow(b, steps=400))-1 for b in b0s])
ax[1].plot(b0s, N, "C0o", ms=3, label="MK steps to strong coupling")
ax[1].plot(b0s, b0s/cAF + (N[-1]-b0s[-1]/cAF), "C3--", label="AF slope $1/0.2575$")
ax[1].set_xlabel("initial Wilson $\\beta$"); ax[1].set_ylabel("$\\log_2(\\xi/a)$"); ax[1].legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig("../figures/mk_flow.pdf")
print("min decrement after first step:", D.min(), " at beta", B[np.argmin(D)])
big = B > 10; print("mean decrement beta>10:", D[big].mean(), " beta>20:", D[B > 20].mean())
fit = np.polyfit(1/B[big], D[big], 1); print("fit delta = %.4f + %.3f/beta" % (fit[1], fit[0]))
print("max decrement:", D.max(), " at beta_eff", B[np.argmax(D)])
print("slope of N vs beta0 (large):", np.polyfit(b0s[b0s > 10], N[b0s > 10], 1)[0], " 1/cAF =", 1/cAF, " 4 =", 4)
