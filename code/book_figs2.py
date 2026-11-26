"""Colour figures for the expanded book (analytic and light numerics). Polson and Sokolov, Nov 2026."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.special import ive
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
out = "../figures/"; rng = np.random.default_rng(2)
# 1 Thorin density heatmap for the toy spectrum
atoms = [(1.0, 1.0), (1.5, 0.6), (1.9, 0.3)]; Ec = np.linspace(2.0005, 25, 1500); rc = 0.4*np.sqrt(Ec-2)*np.exp(-0.3*(Ec-2)); dE = Ec[1]-Ec[0]
al = np.linspace(0.3, 8, 220); tt = np.linspace(0.05, 12, 300)
E_all = np.concatenate([[a for a, _ in atoms], Ec]); w_all = np.concatenate([[w for _, w in atoms], rc*dE])
T = np.array([[np.sum(w_all*np.exp(-a*E_all)*(1-np.cos(t*E_all)))/(np.pi*t) for t in tt] for a in al])
fig, ax = plt.subplots(figsize=(6, 3.4))
im = ax.imshow(T.T, origin="lower", aspect="auto", extent=[al[0], al[-1], tt[0], tt[-1]], norm=LogNorm(vmin=T.max()*1e-4, vmax=T.max()), cmap="magma")
fig.colorbar(im, ax=ax, label="$T_\\alpha$"); ax.set_xlabel("Euclidean time $\\alpha$"); ax.set_ylabel("real time $t=\\sqrt{2z}$")
fig.tight_layout(); fig.savefig(out+"thorin_heatmap.pdf"); plt.close()
# 2 heat kernel on SU(2) as a heatmap in (t, theta)
th = np.linspace(1e-3, np.pi-1e-3, 400); n = np.arange(1, 500)[:, None]; ts = np.geomspace(0.02, 5, 250)
Kh = np.array([((n*np.exp(-t*(n**2-1)/4))*np.sin(n*th[None, :])).sum(0)/np.sin(th)*(2/np.pi)*np.sin(th)**2 for t in ts])
fig, ax = plt.subplots(figsize=(6, 3.3))
im = ax.pcolormesh(th, ts, Kh, shading="auto", cmap="viridis", norm=LogNorm(vmin=1e-3, vmax=Kh.max())); ax.set_yscale("log")
fig.colorbar(im, ax=ax, label="density of $\\theta$"); ax.set_xlabel("rotation angle $\\theta$"); ax.set_ylabel("Brownian time $t$")
fig.tight_layout(); fig.savefig(out+"heat_kernel_heatmap.pdf"); plt.close()
# 3 free field Kallen Lehmann densities for phi^2 at zero momentum, d = 2, 3, 4, with offset
m = 1.0; E = np.linspace(2.0005, 10, 800)
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.9))
for d, c in [(2, "C0"), (3, "C2"), (4, "C3")]:
    rho = (E**2/4-m**2)**((d-3)/2)/E; ax[0].plot(E, rho/rho[200], c, label="$d=%d$" % d)
    a_ = np.linspace(0.3, 6, 120); Cc = [np.trapezoid(rho*np.exp(-a*E), E) for a in a_]
    ax[1].semilogy(a_, Cc/Cc[0], c, label="$d=%d$" % d)
ax[0].set_xlabel("$E$ (units of $m$)"); ax[0].set_ylabel("spectral density (scaled)"); ax[0].set_ylim(0, 3); ax[0].legend(frameon=False, fontsize=8)
ax[1].set_xlabel("$\\tau$"); ax[1].set_ylabel("$C(\\tau)/C(0.3)$"); ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"free_field_kl.pdf"); plt.close()
# 4 Gross Witten Wadia eigenvalue densities
thc = np.linspace(-np.pi, np.pi, 800); fig, ax = plt.subplots(figsize=(5.6, 3.1)); cm = plt.cm.plasma
for k, lam in enumerate([0.5, 1.0, 1.6, 2.0, 3.0, 5.0]):
    if lam >= 2: rho = (1/(2*np.pi))*(1+(2/lam)*np.cos(thc))
    else:
        s2 = np.sin(thc/2)**2; rho = np.where(s2 < lam/2, (2/(np.pi*lam))*np.cos(thc/2)*np.sqrt(np.clip(lam/2-s2, 0, None)), 0)
    ax.plot(thc, rho, color=cm(k/6), label="$\\lambda=%g$" % lam)
ax.set_xlabel("eigenvalue angle"); ax.set_ylabel("density"); ax.legend(frameon=False, fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(out+"gww.pdf"); plt.close()
# 5 Polya random walk: fraction of walkers returned by time n, d = 1, 2, 3
fig, ax = plt.subplots(figsize=(5.4, 3.1)); nmax = 4000; W = 4000
for d, c in [(1, "C0"), (2, "C1"), (3, "C2")]:
    pos = np.zeros((W, d), int); ret = np.zeros(W, bool); frac = []
    for s in range(nmax):
        ax_ = rng.integers(d, size=W); sg = rng.choice([-1, 1], size=W); pos[np.arange(W), ax_] += sg
        ret |= np.all(pos == 0, axis=1); frac.append(ret.mean())
    ax.semilogx(np.arange(1, nmax+1), frac, c, label="$d=%d$" % d)
ax.axhline(0.3405, color="C2", ls=":", lw=1, label="P\\'olya constant $0.3405$")
ax.set_xlabel("steps $n$"); ax.set_ylabel("fraction returned to origin"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"polya_walk.pdf"); plt.close()
# 6 Ehrenfest urn: Kac's lemma, mean return time = 1/pi(k)
from math import comb
Nb = 20; pi = np.array([comb(Nb, k) for k in range(Nb+1)])/2**Nb
x = Nb//2; visits = {}; last = {}; rts = {k: [] for k in range(Nb+1)}
for s in range(2_000_000):
    x = x-1 if rng.random() < x/Nb else x+1
    if x in last: rts[x].append(s-last[x])
    last[x] = s
ks = np.arange(Nb+1); mrt = np.array([np.mean(rts[k]) if len(rts[k]) > 5 else np.nan for k in ks])
fig, ax = plt.subplots(figsize=(5.4, 3.1))
ax.semilogy(ks, 1/pi, "k-", lw=2, alpha=0.4, label="Kac: $1/\\pi(k)$")
ax.semilogy(ks, mrt, "C3o", ms=4, label="simulated mean return time")
ax.set_xlabel("balls in urn A, $k$"); ax.set_ylabel("mean recurrence time"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"ehrenfest_kac.pdf"); plt.close()
# 7 dynamic versus static gap: Ising ring, heat bath Glauber generator spectral gap
nn = 10; states = np.array([[1 if (s >> i) & 1 else -1 for i in range(nn)] for s in range(2**nn)])
bet = np.linspace(0.05, 1.6, 24); dyn = []
for b in bet:
    Q = np.zeros((2**nn, 2**nn))
    for s in range(2**nn):
        for i in range(nn):
            t_ = s ^ (1 << i); sp = states[s]; h = sp[(i-1) % nn]+sp[(i+1) % nn]
            rate = 1/(1+np.exp(2*b*sp[i]*h))/nn      # heat bath flip rate
            Q[s, t_] = rate
        Q[s, s] = -Q[s].sum()
    p = np.exp(b*np.sum(states*np.roll(states, 1, axis=1), axis=1)); p /= p.sum()
    D = np.sqrt(p); Ssym = (D[:, None]*Q)/D[None, :]; ev = np.sort(np.linalg.eigvalsh((Ssym+Ssym.T)/2))[::-1]
    dyn.append(-ev[1])
stat = -np.log(np.tanh(bet))
fig, ax = plt.subplots(figsize=(5.4, 3.1))
ax.semilogy(bet, dyn, "C0o-", ms=3, label="Glauber spectral gap (dynamic)")
ax.semilogy(bet, stat, "C3-", label="transfer matrix gap $-\\log\\tanh\\beta$ (static)")
ax.set_xlabel("$\\beta$"); ax.legend(frameon=False, fontsize=8); ax.set_ylabel("gap")
fig.tight_layout(); fig.savefig(out+"dynamic_static_gap.pdf"); plt.close()
np.save("/tmp/dynstat.npy", np.vstack([bet, dyn, stat]))
# 8 Gaussian ledger: flat versus curved precision flow (heat kernel coordinate)
fig, ax = plt.subplots(figsize=(5.4, 3.1)); cm = plt.cm.viridis
for k, b0 in enumerate([2, 4, 6, 8, 10]):
    bb = [b0]; 
    for _ in range(40): bb.append(bb[-1]-0.25 if bb[-1] > 0.6 else bb[-1]*0.5)
    ax.plot(range(41), bb, color=cm(k/5), lw=1.8, label="curved, $\\beta_0=%d$" % b0)
    ax.plot(range(41), [b0]*41, color=cm(k/5), ls=":", lw=1)
ax.set_xlabel("renormalization step"); ax.set_ylabel("coupling $\\beta$"); ax.legend(frameon=False, fontsize=7, ncol=2)
ax.text(30, 9.4, "flat space: marginal", fontsize=8, color="gray")
fig.tight_layout(); fig.savefig(out+"gaussian_ledger.pdf"); plt.close()
print("ok")
