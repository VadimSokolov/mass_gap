"""Yakowitz, Krimmel and Szidarovszky weighted Monte Carlo: trapezoid rule on sorted uniform points,
mean square error O(n^-4) for C^2 integrands (Masry and Cambanis), against plain Monte Carlo O(n^-1).
Integrand: <S>_b for 2D SU(2) on the 3x3 torus, exact; integral -int_0^2 <S>_b db = log Z(2).
Polson and Sokolov, Nov 2026."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import vl_core as v
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
bg = np.linspace(0, 2, 801); ms = np.array([float(v.P) if b == 0 else v.meanS_exact(b, h=min(1e-5, b/4)) for b in bg]); f = CubicSpline(bg, ms)
exact = -float(f.integrate(0, 2)); print("check vs exact log Z(2):", exact, v.logZ_exact(2.0))
rng = np.random.default_rng(1); ns = [4, 8, 16, 32, 64, 128, 256]; rm, ry = [], []
for n in ns:
    em, ey = [], []
    for _ in range(3000):
        u = np.sort(rng.random(n))*2
        em.append(-2*np.mean(f(u)) - exact)
        x = np.concatenate([[0], u, [2]]); ey.append(-np.trapezoid(f(x), x) - exact)
    rm.append(np.sqrt(np.mean(np.square(em)))); ry.append(np.sqrt(np.mean(np.square(ey))))
    print("n=%4d  RMSE plain MC %.2e   Yakowitz trapezoid %.2e" % (n, rm[-1], ry[-1]))
sm = np.polyfit(np.log(ns), np.log(rm), 1)[0]; sy = np.polyfit(np.log(ns[2:]), np.log(ry[2:]), 1)[0]
print("slopes: plain %.2f (theory -0.5), Yakowitz %.2f (theory -2)" % (sm, sy))
fig, ax = plt.subplots(figsize=(5.2, 3.2))
ax.loglog(ns, rm, "C0o-", label="plain Monte Carlo, slope %.2f" % sm)
ax.loglog(ns, ry, "C3s-", label="Yakowitz trapezoid, slope %.2f" % sy)
ax.set_xlabel("number of random rungs $n$"); ax.set_ylabel("RMSE of $\\log Z(2)$"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("../figures/yakowitz_ti.pdf"); print("ok")
