# Three dimensional ladder: beta a sqrt(sigma) as a function of log2(beta), against lattice and Karabali Nair.
import numpy as np, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "code"); from mk_sigma import sigma
lb = np.linspace(1, 5, 81); v = np.array([2**x*np.sqrt(sigma(2**x, d=3, steps=80)) for x in lb])
print("min %.3f max %.3f mean over one period %.3f" % (v[40:].min(), v[40:].max(), v[60:80].mean()))
fig, ax = plt.subplots(figsize=(5.6, 3.3)); sc = ax.scatter(lb, v, c=lb % 1, cmap="twilight", s=14); ax.plot(lb, v, color="#7f8c8d", lw=0.7)
ax.axhline(1.341, color="#c0392b", label="lattice, 1.341"); ax.axhline(4*np.sqrt(3/(8*np.pi)), color="#e67e22", ls="--", label="Karabali, Kim and Nair, 1.382")
ax.set_ylim(0.5, 1.6); ax.set_xlabel(r"$\log_2\beta$"); ax.set_ylabel(r"$\beta\,a\sqrt{\sigma}$"); ax.legend(frameon=False, fontsize=7); fig.colorbar(sc, label="phase of $\\log_2\\beta$")
fig.tight_layout(); fig.savefig("figures/mk_3d.pdf")
