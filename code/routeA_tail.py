"""Route A: the Yang Mills GGC X_alpha for the toy spectrum, simulated as a gamma convolution.
Atoms of the Thorin measure at rate E with mass w give Gamma(shape w, rate E). Tail slope = -Delta."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(3); alpha = 0.5
atoms = [(1.0, 1.0), (1.5, 0.6), (1.9, 0.3)]
E = np.linspace(2.0005, 30, 3000); dE = E[1]-E[0]
rho = 0.4*np.sqrt(E-2)*np.exp(-0.3*(E-2))
rates = np.array([a for a, _ in atoms] + list(E))
mass = np.array([w for _, w in atoms] + list(rho*dE))*np.exp(-alpha*rates)
N = 2_000_000
X = np.zeros(N)
for r, m in zip(rates, mass):
    if m > 1e-12: X += rng.gamma(m, 1/r, N)
xs = np.linspace(0, 14, 200); tail = np.array([(X > x).mean() for x in xs])
fig, ax = plt.subplots(figsize=(5.2, 3.0))
ax.semilogy(xs, tail, "C0", label="$\\Pr(X_\\alpha>x)$, simulated")
ax.semilogy(xs, tail[60]*np.exp(-(xs-xs[60])), "k:", label="slope $-\\Delta=-1$")
ax.set_xlabel("$x$"); ax.legend(frameon=False, fontsize=8); ax.set_ylim(1e-6, 1.2)
fig.tight_layout(); fig.savefig("../figures/routeA_tail.pdf")
m = (xs > 4) & (tail > 1e-5)
print("fitted tail slope:", np.polyfit(xs[m], np.log(tail[m]), 1)[0])
