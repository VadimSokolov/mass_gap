import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
K = np.linspace(0.05, 1.2, 400)
Ks = -0.5*np.log(np.tanh(K))          # dual coupling
gap2d = 2*np.abs(K-Ks); Kc = 0.5*np.log(1+np.sqrt(2))
gap1d = -np.log(np.tanh(K))
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.9))
ax[0].plot(K, gap1d, "C0", label="Ising chain $-\\log\\tanh K$")
ax[0].plot(K, 2*np.exp(-2*K), "C0:", label="$2e^{-2K}$")
ax[0].set_ylim(0, 4); ax[0].set_xlabel("$K$"); ax[0].set_ylabel("gap (lattice units)"); ax[0].legend(frameon=False, fontsize=7.5)
ax[1].plot(K, gap2d, "C3", label="2D Ising $2|K-K^*|$")
ax[1].axvline(Kc, color="k", ls=":", label="$K_c=\\frac{1}{2}\\log(1+\\sqrt{2})$")
ax[1].set_ylim(0, 3); ax[1].set_xlabel("$K$"); ax[1].legend(frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig("../figures/examples_gaps.pdf")
# massless example: mu(dE) = E^{kappa-1} dE, Route B local rate tends to 0
al = np.linspace(0.5, 40, 200)
fig, ax = plt.subplots(figsize=(5.2, 2.8))
for kap, c in [(0.5, "C0"), (1.0, "C1"), (2.0, "C2")]:
    E = np.linspace(1e-4, 60, 60000); t = 1.3
    T = np.array([np.sum(E**(kap-1)*np.exp(-a*E)*(1-np.cos(t*E)))*(E[1]-E[0]) for a in al])
    ax.plot(al[1:], -np.diff(np.log(T))/np.diff(al), c, label="$\\kappa=%g$" % kap)
ax.set_xlabel("$\\alpha$"); ax.set_ylabel("local Thorin rate"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("../figures/massless_rates.pdf")
print(Kc, "ok")
