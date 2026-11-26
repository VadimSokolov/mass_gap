# Kadanoff's variational bound: the Migdal Kadanoff ladder gives a rigorous upper bound on the pressure
#   P(beta) = lim (1/#plaquettes) log int prod_p exp(beta cos theta_p) dHaar,   P'(beta) = W11(beta),
#   P(beta) <= P_MK(beta) = log z0 + sum_m b^(-d m) (1/k) log int g_m^k dHaar,   k = b^(d-2),
# with g_m the normalised plaquette law at rung m. By convexity, W11(beta) <= inf_{beta' > beta} P_MK(beta')/(beta' - beta).
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ive
from scipy.interpolate import PchipInterpolator
M, NMAX = 3000, 300; x, w = np.polynomial.legendre.leggauss(M); th = (x+1)*np.pi/2; w = w*np.pi/2
haar = (2/np.pi)*np.sin(th)**2*w; n = np.arange(1, NMAX+1); CHI = np.sin(np.outer(n, th))/np.sin(th)
def P_mk(beta, d=4, b=2, steps=30):
    k = b**(d-2); u = np.clip(ive(n, beta)/ive(1, beta), 1e-300, None); P = np.log(2*ive(1, beta)/beta)+beta; fac = 1.0
    for m in range(steps):
        g = np.clip((n*u)@CHI, 1e-300, None); c = (g**k*haar).sum(); P += fac*np.log(c)/k; fac /= b**d
        cc = (CHI@(g**k*haar))/n; u = np.clip(cc/cc[0], 1e-300, None)**(b*b)
        if fac < 1e-14: break
    return P
# measured plaquette of four dimensional SU(2), 8^4 lattice (this book), plus the strong coupling value beta/4 near zero
bm = np.array([0, 0.2, 0.4, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 4.0, 6.0, 10.0])
wm = np.array([0, 0.05, 0.0998, 0.2901, 0.3380, 0.3875, 0.4411, 0.5014, 0.5700, 0.6527, 0.7237, 0.7996, 0.8700, 0.9233])
W = PchipInterpolator(bm, wm); P_true = lambda b: W.integrate(0, b)
grid = np.linspace(0.05, 14, 280); Pm = np.array([P_mk(b) for b in grid])
print(" beta   P_true (integrated plaquette)   P_MK (rigorous upper bound)   P_2D = log z0 + beta")
for b in [1.0, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0]: print(" %4.1f   %.4f   %.4f   %.4f" % (b, P_true(b), P_mk(b), np.log(2*ive(1, b)/b)+b))
print(" beta   measured W11   ladder bound inf P_MK(b')/(b'-b)   heat bath bound u(6 beta)")
rows = []
for b in [1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0]:
    sel = grid > b+0.05; bound = (Pm[sel]/(grid[sel]-b)).min(); rows.append((b, float(W(b)), bound, ive(2, 6*b)/ive(1, 6*b))); print(" %4.1f   %.4f   %.4f   %.4f" % rows[-1])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1)); bb = np.linspace(0.05, 10, 200)
ax[0].plot(bb, [P_true(b) for b in bb], color="#c0392b", label="true, integrated plaquette"); ax[0].plot(grid, Pm, color="#1f4e79", label="ladder, upper bound")
ax[0].plot(bb, np.log(2*ive(1, bb)/bb)+bb, color="#16a085", ls=":", label="two dimensions, exact"); ax[0].set_xlim(0, 10); ax[0].set_ylim(0, 9)
ax[0].set_xlabel(r"$\beta$"); ax[0].set_ylabel("pressure per plaquette"); ax[0].legend(frameon=False, fontsize=7)
ax[1].plot(bb, np.interp(bb, grid, Pm)-np.array([P_true(b) for b in bb]), color="#8e44ad"); ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"$\beta$"); ax[1].set_ylabel(r"$P_{\mathrm{MK}}-P$"); fig.tight_layout(); fig.savefig("figures/kadanoff_bound.pdf")
