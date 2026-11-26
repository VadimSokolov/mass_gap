# Further colour figures: spectrum and gap, Wilson laws, Brownian motion on S^3, reflection positivity, 4D histories.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.special import iv
F = "figures/"; rng = np.random.default_rng(11)
# 1. spectrum with a gap and the correlators it produces
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
E = np.linspace(0, 6, 600); atoms = [(1.0, 0.6), (1.6, 0.25)]
cont = np.where(E > 2.0, 0.25*np.sqrt(np.clip(E-2, 0, None))*np.exp(-0.8*(E-2)), 0)
ax[0].fill_between(E, cont, color="#5dade2", alpha=0.7, label="continuum")
for e, w in atoms: ax[0].vlines(e, 0, w, color="#c0392b", lw=3)
ax[0].axvspan(0, 1.0, color="#f5b041", alpha=0.25); ax[0].text(0.2, 0.5, r"gap $\Delta$", fontsize=9)
ax[0].set_xlabel("energy $E$"); ax[0].set_ylabel("spectral weight"); ax[0].set_title("spectral measure", fontsize=9)
tau = np.linspace(0, 8, 300); col = cm.viridis(np.linspace(0.1, 0.9, 4))
for D, c in zip([0.5, 1.0, 1.5, 2.0], col):
    C = sum(w*np.exp(-(e-1+D)*tau) for e, w in atoms) + np.trapezoid(cont[None, :]*np.exp(-np.outer(tau, E+D-1)), E, axis=1)
    ax[1].semilogy(tau, C/C[0], color=c, label=fr"$\Delta={D}$")
ax[1].set_xlabel(r"Euclidean time $\tau$"); ax[1].set_ylabel("correlation"); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title("decay at rate $\\Delta$", fontsize=9)
fig.tight_layout(); fig.savefig(F+"spectrum_gap.pdf"); plt.close()
# 2. Wilson plaquette laws on SU(2) against Haar
th = np.linspace(0, np.pi, 400); fig, ax = plt.subplots(figsize=(5.4, 3.4)); bs = [0.5, 1, 2, 4, 8, 16]
for b, c in zip(bs, cm.plasma(np.linspace(0.05, 0.85, len(bs)))):
    d = (2/np.pi)*np.sin(th)**2*np.exp(b*np.cos(th))/ (2*iv(1, b)/b); ax.plot(th, d, color=c, label=fr"$\beta={b}$")
ax.plot(th, (2/np.pi)*np.sin(th)**2, "k--", lw=1, label="Haar")
ax.set_xlabel(r"rotation angle $\theta$ of $U_p$"); ax.set_ylabel("density"); ax.legend(frameon=False, fontsize=7, ncol=2)
ax.set_title("Wilson plaquette laws on SU(2)", fontsize=9); fig.tight_layout(); fig.savefig(F+"wilson_laws.pdf"); plt.close()
# 3. Brownian motion on S^3: angle from identity for several paths
def qm(a, b):
    a0,a1,a2,a3 = a; b0,b1,b2,b3 = b
    return np.array([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0])
dt = 0.002; nT = 3000; fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [2.2, 1]})
finals = []
for j, c in enumerate(cm.turbo(np.linspace(0.05, 0.95, 8))):
    g = np.array([1., 0, 0, 0]); ang = []
    for i in range(nT):
        v = rng.standard_normal(3)*np.sqrt(dt); h = np.concatenate([[np.cos(np.linalg.norm(v))], np.sinc(np.linalg.norm(v)/np.pi)*v])
        g = qm(g, h); g /= np.linalg.norm(g); ang.append(np.arccos(np.clip(g[0], -1, 1)))
    ax[0].plot(np.arange(nT)*dt, ang, color=c, lw=0.8)
for k in range(3000):
    g = np.array([1., 0, 0, 0])
    v = rng.standard_normal((300, 3))*np.sqrt(0.02)
    for vv in v:
        n = np.linalg.norm(vv); g = qm(g, np.concatenate([[np.cos(n)], np.sinc(n/np.pi)*vv])); g /= np.linalg.norm(g)
    finals.append(np.arccos(np.clip(g[0], -1, 1)))
ax[0].set_xlabel("Brownian time"); ax[0].set_ylabel(r"angle $\theta$ from identity"); ax[0].set_title("paths on SU(2)", fontsize=9)
ax[1].hist(finals, bins=30, density=True, color="#76b7b2", orientation="horizontal")
ax[1].plot((2/np.pi)*np.sin(th)**2, th, "k--", lw=1); ax[1].set_xlabel(r"density at time 6 ($K_t$ time 12)"); ax[1].set_title("approach to Haar", fontsize=9)
fig.tight_layout(); fig.savefig(F+"s3_brownian.pdf"); plt.close()
# 4. reflection positivity: Gram matrix of a correlator and of a non positive function
tau = np.arange(0, 12)*0.4
C = 0.6*np.exp(-1.0*tau) + 0.3*np.exp(-1.7*tau) + 0.1*np.exp(-2.6*tau)
B = np.exp(-tau)*np.cos(2.2*tau)
G1 = C[:6, None]*0; G1 = np.array([[C[i+j] for j in range(6)] for i in range(6)]); G2 = np.array([[B[i+j] for j in range(6)] for i in range(6)])
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.7))
for a, G, t in [(ax[0], G1, "reflection positive"), (ax[1], G2, r"$e^{-\tau}\cos 2.2\tau$")]:
    im = a.imshow(G, cmap="RdBu_r", vmin=-1, vmax=1); a.set_title(t, fontsize=8); a.set_xticks([]); a.set_yticks([])
fig.colorbar(im, ax=ax[:2], fraction=0.03)
ax[2].semilogy(np.arange(1, 7), np.sort(np.abs(np.linalg.eigvalsh(G1)))[::-1], "o-", color="#1f4e79", label="positive")
e2 = np.sort(np.linalg.eigvalsh(G2))[::-1]; ax[2].plot(np.arange(1, 7), np.abs(e2), "s-", color="#c0392b", label="abs, not positive")
ax[2].scatter(np.arange(1, 7)[e2 < -1e-12], np.abs(e2[e2 < -1e-12]), s=80, facecolors="none", edgecolors="k")
ax[2].set_title("Hankel eigenvalues", fontsize=8); ax[2].legend(frameon=False, fontsize=6); ax[2].set_xlabel("index")
fig.savefig(F+"rp_hankel.pdf", bbox_inches="tight"); plt.close()
# 5. 4D plaquette histories
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0)); bs = ["2.2", "2.3", "2.4", "2.5", "2.635"]
for b, c in zip(bs, cm.viridis(np.linspace(0.05, 0.9, 5))):
    P = np.load(f"/tmp/plaq_hist_{b}.npy"); ax[0].plot(P, color=c, lw=0.6, label=fr"$\beta={b}$")
    x = P-P.mean(); r = [1]+[np.dot(x[:-k], x[k:])/np.dot(x, x) for k in range(1, 30)]; ax[1].plot(r, "o-", ms=2, color=c)
ax[0].set_xlabel("sweep"); ax[0].set_ylabel("mean plaquette"); ax[0].legend(frameon=False, fontsize=6, ncol=2)
ax[1].axhline(0, color="k", lw=0.5); ax[1].set_xlabel("lag"); ax[1].set_ylabel("autocorrelation")
fig.tight_layout(); fig.savefig(F+"plaquette_history.pdf"); plt.close()
print("ok", np.linalg.eigvalsh(G1).min(), np.linalg.eigvalsh(G2).min())
