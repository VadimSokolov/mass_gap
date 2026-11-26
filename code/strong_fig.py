import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ive
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
b = np.linspace(0.05, 4, 300); u = ive(2, b)/ive(1, b)
fig, ax = plt.subplots(figsize=(5.4, 3.1))
ax.semilogy(b, -4*np.log(u), "C0", label="leading strong coupling $-4\\log u_2(\\beta)$")
bb = np.linspace(1.5, 4, 100)
for c, ls in [(30, "--"), (100, ":")]:
    ax.semilogy(bb, c*np.exp(-3*np.pi**2*bb/11), "C3"+ls if ls == "--" else "C3:", label="$%d\\,e^{-3\\pi^2\\beta/11}$ (scaling, illustrative)" % c)
ax.set_xlabel("$\\beta$ (SU(2) Wilson)"); ax.set_ylabel("lattice mass $m(\\beta)$"); ax.legend(frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig("../figures/strong_vs_scaling.pdf"); print("ok")
for x in [0.1, 0.5, 1.0, 1.5, 2.0]: print(x, ive(2,x)/ive(1,x), -4*np.log(ive(2,x)/ive(1,x)))
