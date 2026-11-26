"""Figure for open problem 1. Polson and Sokolov, Nov 2026."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ive
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
x = np.logspace(-2.5, 1.2, 400); beta = 1.5/x
ratio = -beta*np.log(ive(2, beta)/ive(1, beta))/1.5
fig, ax = plt.subplots(1, 2, figsize=(10, 3.7))
ax[0].semilogx(x, ratio, "C0", label="$\\Delta(a)/\\Delta_{\\mathrm{cont}}$ (exact)")
ax[0].semilogx(x, np.arcsinh(x)/x, "k--", label="bound $\\operatorname{arcsinh}(x)/x$")
ax[0].axhline(1, color="gray", lw=0.6)
ax[0].set_xlabel("$x = 3g^2a^2/8$"); ax[0].legend(frameon=False); ax[0].set_ylim(0, 1.25)
# Thorin local rates in physical units, operator O = chi_2+chi_3+chi_4, g = l = 1
g, l = 1.0, 3.0; dc = 3*g*g*l/8
for L, col in [(1, "C3"), (2, "C1"), (4, "C2"), (64, "C0")]:
    a = l/L; b = 4/(g*g*a*a)
    E = np.array([-(l/a**2)*np.log(ive(n, b)/ive(1, b)) for n in (2, 3, 4)])
    t = 2.0; al = np.linspace(0.2, 25, 300)
    T = np.array([np.sum(np.exp(-A*E)*2*np.sin(t*E/2)**2) for A in al])
    ax[1].plot(al[1:], -np.diff(np.log(T))/np.diff(al)/dc, col, label="L = %d" % L)
ax[1].axhline(1, color="k", ls=":", label="$\\Delta_{\\mathrm{cont}}$")
ax[1].axhline(np.arcsinh(3*g*g*l*l/8)/(3*g*g*l*l/8), color="k", ls="--", lw=0.8, label="uniform $\\Delta_0$")
ax[1].set_ylim(0.2, 2.5); ax[1].set_xlabel("$\\alpha$ (physical Euclidean time)")
ax[1].set_ylabel("$-\\partial_\\alpha\\log T_\\alpha / \\Delta_{\\mathrm{cont}}$"); ax[1].legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig("../figures/op1_2d_uniform.pdf"); print("ok")
