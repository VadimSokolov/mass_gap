# Empirical sandwich: the law of the 2 x 2 loop holonomy against two caricatures built from the measured plaquette law:
# independent plaquettes (decimation only) and Migdal Kadanoff (bond moving with power 4, then decimation).
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
bins = np.linspace(0, np.pi, 61); th = 0.5*(bins[1:]+bins[:-1]); dth = bins[1]-bins[0]; n = np.arange(1, 40)
CHI = np.sin(np.outer(n, th))/np.sin(th); haar = (2/np.pi)*np.sin(th)**2*dth
def coeff(dens): c = (CHI@(dens*dth))/n; return c/c[0]            # dens: density of the angle on (0, pi)
def dens_from(u): return np.clip((n*u)@CHI, 1e-12, None)*(2/np.pi)*np.sin(th)**2
def monotone(r, wts):                                              # weighted share of bins where the ratio increases
    ok = np.isfinite(r[1:]) & np.isfinite(r[:-1]); d = np.diff(r)[ok]; w = np.minimum(wts[1:], wts[:-1])[ok]; return (w*(d > 0)).sum()/w.sum()
fig, ax = plt.subplots(2, 3, figsize=(7.4, 4.6), sharex=True); out = []
for a, b in zip(ax.ravel(), ["2.2", "2.5", "3.0", "4.0", "6.0", "10.0"]):
    h1, h2 = np.load(f"/tmp/loops_{b}.npy"); p1 = h1/h1.sum()/dth; p2 = h2/h2.sum()/dth
    u = coeff(p1); dec = dens_from(u**4)
    f = np.clip((n*u)@CHI, 1e-12, None); g = f**4; c = (CHI@(g*haar))/n; mk = dens_from((c/c[0])**4)
    a.plot(th, dec, color="#16a085", label="independent"); a.plot(th, mk, color="#1f4e79", label="Migdal Kadanoff"); a.plot(th, p2, ".", ms=3, color="#c0392b", label="measured")
    a.set_title(fr"$\beta={b}$", fontsize=8); m1 = monotone(dec/np.where(p2 > 0, p2, np.nan), p2); m2 = monotone(p2/mk, p2)
    out.append((float(b), coeff(dec)[1], coeff(p2)[1], coeff(mk)[1], m1, m2))
    print("beta %5s  u2: independent %.4f  measured %.4f  MK %.4f   ratio monotone shares %.2f %.2f" % ((b,)+out[-1][1:]))
ax[0, 0].legend(frameon=False, fontsize=6); [x.set_xlabel(r"$\theta$") for x in ax[1]]; [x.set_ylabel("density") for x in ax[:, 0]]
fig.tight_layout(); fig.savefig("figures/sandwich.pdf")
