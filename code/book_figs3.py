"""Simulation figures: Jerrum Sinclair style product estimator, Swendsen Wang vs Metropolis with Madras Sokal
windowing, topological freezing in 2D U(1), Lee Yang and Fisher zeros. Polson and Sokolov, Nov 2026."""
import numpy as np, matplotlib, itertools
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
out = "../figures/"; rng = np.random.default_rng(7)
def tau_int(x, c=6.0):
    """Madras Sokal automatic windowing estimate of the integrated autocorrelation time"""
    x = np.asarray(x, float)-np.mean(x); n = len(x); f = np.fft.rfft(x, 2*n); acf = np.fft.irfft(f*np.conj(f))[:n]; acf /= acf[0]
    tau = 0.5
    for W in range(1, n//2):
        tau += acf[W]
        if W >= c*tau: break
    return max(tau, 0.5)
# Ising helpers on an L x L torus
def energy(s): return -np.sum(s*np.roll(s, 1, 0)) - np.sum(s*np.roll(s, 1, 1))
def metro_sweep(s, b):
    Lx = s.shape[0]
    for par in (0, 1):
        mask = (np.add.outer(np.arange(Lx), np.arange(Lx)) % 2) == par
        h = np.roll(s, 1, 0)+np.roll(s, -1, 0)+np.roll(s, 1, 1)+np.roll(s, -1, 1)
        dE = 2*s*h; acc = (rng.random(s.shape) < np.exp(-b*dE)) & mask; s[acc] *= -1
    return s
def sw_step(s, b):
    Lx = s.shape[0]; N = Lx*Lx; idx = np.arange(N).reshape(Lx, Lx); p = 1-np.exp(-2*b)
    rows, cols = [], []
    for ax in (0, 1):
        nb = np.roll(idx, -1, ax); same = s == np.roll(s, -1, ax); on = same & (rng.random(s.shape) < p)
        rows.append(idx[on]); cols.append(nb[on])
    r = np.concatenate(rows); c = np.concatenate(cols)
    G = coo_matrix((np.ones(len(r)), (r, c)), shape=(N, N)); nc, lab = connected_components(G, directed=False)
    flip = rng.random(nc) < 0.5; s = s.reshape(-1); s[flip[lab]] *= -1; return s.reshape(Lx, Lx)
# 1 exact 4x4 Ising by enumeration: energy and magnetisation histogram
L4 = 4; confs = np.array(list(itertools.product([-1, 1], repeat=16))).reshape(-1, 4, 4)
Es = np.array([energy(c) for c in confs]); Ms = confs.reshape(len(confs), -1).sum(1)
logZ = lambda b: np.log(np.sum(np.exp(-b*Es)))
# 2 Jerrum Sinclair style product estimator Z(b)/Z(0) = prod E_{b_i}[exp(-(b_{i+1}-b_i) E)], SW sampler
bfin = 0.6; K = 12; sched = np.linspace(0, bfin, K+1); exact = logZ(bfin) - 16*np.log(2)
errs = {}
for nsamp in [20, 80, 320]:
    e = []
    for rep in range(60):
        tot = 0.0
        for i in range(K):
            s = rng.choice([-1, 1], size=(4, 4)); vals = []
            for _ in range(10): s = sw_step(s, sched[i])
            for _ in range(nsamp): s = sw_step(s, sched[i]); vals.append(np.exp(-(sched[i+1]-sched[i])*energy(s)))
            tot += np.log(np.mean(vals))
        e.append(tot-exact)
    errs[nsamp] = np.array(e); print("product estimator n=%d: rmse %.4f" % (nsamp, np.sqrt(np.mean(np.square(e)))))
fig, ax = plt.subplots(figsize=(5.4, 3.1)); cols = ["C0", "C1", "C3"]
for c, (k, e) in zip(cols, errs.items()): ax.hist(e, bins=20, alpha=0.6, color=c, label="%d samples per rung" % k)
ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("error in $\\log Z(0.6)$, $4\\times4$ Ising"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"product_estimator.pdf"); plt.close()
# 3 autocorrelation: Metropolis vs Swendsen Wang on 24 x 24, energy observable, Madras Sokal windowing
Lb = 24; bc = 0.5*np.log(1+np.sqrt(2)); bets = [0.30, 0.38, 0.42, bc, 0.46, 0.52]; tm, ts = [], []
for b in bets:
    s = rng.choice([-1, 1], size=(Lb, Lb)); xs = []
    for _ in range(300): s = metro_sweep(s, b)
    for _ in range(4000): s = metro_sweep(s, b); xs.append(abs(s.sum()))
    tm.append(tau_int(xs))
    s = rng.choice([-1, 1], size=(Lb, Lb)); xs = []
    for _ in range(100): s = sw_step(s, b)
    for _ in range(2000): s = sw_step(s, b); xs.append(abs(s.sum()))
    ts.append(tau_int(xs)); print("beta %.3f tau Metropolis %.1f  SW %.1f" % (b, tm[-1], ts[-1]))
fig, ax = plt.subplots(figsize=(5.4, 3.1))
ax.semilogy(bets, tm, "C0o-", label="Metropolis (sweeps)"); ax.semilogy(bets, ts, "C3s-", label="Swendsen Wang (steps)")
ax.axvline(bc, color="k", ls=":", lw=1, label="$\\beta_c$"); ax.set_xlabel("$\\beta$"); ax.set_ylabel("$\\tau_{\\mathrm{int}}$ of $|M|$")
ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(out+"sw_vs_metropolis.pdf"); plt.close()
# 4 topological freezing in 2D U(1): theta links, Metropolis, topological charge Q = sum wrapped plaquette angle / 2pi
Lu = 20
def plaq_angles(t):
    return t[:, :, 0]+np.roll(t[:, :, 1], -1, 0)-np.roll(t[:, :, 0], -1, 1)-t[:, :, 1]
def wrap(a): return (a+np.pi) % (2*np.pi)-np.pi
def u1_sweep(t, b, eps):
    for mu in (0, 1):
        for par in (0, 1):
            mask = (np.add.outer(np.arange(Lu), np.arange(Lu)) % 2) == par
            new = t.copy(); new[:, :, mu] = np.where(mask, t[:, :, mu]+eps*rng.uniform(-1, 1, (Lu, Lu)), t[:, :, mu])
            S_old = -b*np.cos(plaq_angles(t)); S_new = -b*np.cos(plaq_angles(new))
            # local action change for each link: plaquettes containing the link
            if mu == 0: dS = (S_new-S_old)+np.roll(S_new-S_old, 1, 1)
            else: dS = (S_new-S_old)+np.roll(S_new-S_old, 1, 0)
            acc = mask & (rng.random((Lu, Lu)) < np.exp(-dS)); t[:, :, mu] = np.where(acc, new[:, :, mu], t[:, :, mu])
    return t
hist = {}; taus = {}
for b in [1.0, 3.0, 6.0, 10.0]:
    t = rng.uniform(-np.pi, np.pi, (Lu, Lu, 2)); eps = 1.5/np.sqrt(b); Qs = []
    for _ in range(200): t = u1_sweep(t, b, eps)
    for _ in range(3000): t = u1_sweep(t, b, eps); Qs.append(np.round(np.sum(wrap(plaq_angles(t)))/(2*np.pi)))
    hist[b] = np.array(Qs); taus[b] = tau_int(Qs); print("U(1) beta %.1f tau_Q %.1f" % (b, taus[b]))
fig, ax = plt.subplots(2, 1, figsize=(6, 4.2), sharex=True); cm = plt.cm.plasma
for k, (b, Q) in enumerate(hist.items()):
    ax[0].plot(Q+0.0, color=cm(k/4), lw=0.7, label="$\\beta=%g$, $\\tau_Q\\approx%.0f$" % (b, taus[b]))
ax[0].set_ylabel("topological charge $Q$"); ax[0].legend(frameon=False, fontsize=7, ncol=2)
for k, (b, Q) in enumerate(hist.items()):
    x = Q-Q.mean(); ac = np.correlate(x, x, "full")[len(x)-1:]; ac /= ac[0]; ax[1].plot(ac[:400], color=cm(k/4))
ax[1].set_xlabel("Monte Carlo sweeps"); ax[1].set_ylabel("autocorrelation of $Q$")
fig.tight_layout(); fig.savefig(out+"u1_topology.pdf"); plt.close()
# 5 Lee Yang zeros (field) and Fisher zeros (temperature) for 4x4 Ising
fig, ax = plt.subplots(1, 2, figsize=(6.4, 3.0)); cm = plt.cm.viridis
for k, b in enumerate([0.2, 0.4, 0.6]):
    # Z = sum_M c_M z^{(16-M)/2}, z = exp(-2h)
    coef = np.zeros(17)
    for E_, M_ in zip(Es, Ms): coef[(16-M_)//2] += np.exp(-b*E_)
    r = np.roots(coef[::-1]); ax[0].plot(r.real, r.imag, "o", ms=3.5, color=cm(k/3), label="$\\beta=%g$" % b)
tc = np.linspace(0, 2*np.pi, 300); ax[0].plot(np.cos(tc), np.sin(tc), "k:", lw=0.6); ax[0].set_aspect("equal")
ax[0].set_title(r"Lee Yang zeros in $e^{-2\eta}$", fontsize=9); ax[0].legend(frameon=False, fontsize=7)
# Fisher: Z = sum_E g(E) x^{(E+32)/4}, x = exp(-4 beta)? energies step 4; use y = exp(-2 beta) with E in steps of 4
coefE = np.zeros(33)
for E_ in Es: coefE[(E_+32)//2] += 1
r = np.roots(coefE[::-1]); r = r[np.abs(r) > 1e-8]
ax[1].plot(r.real, r.imag, "o", ms=3.5, color="C3"); ax[1].axhline(0, color="k", lw=0.4); ax[1].axvline(0, color="k", lw=0.4)
ax[1].set_title("Fisher zeros in $e^{-2\\beta}$", fontsize=9); ax[1].set_aspect("equal")
fig.tight_layout(); fig.savefig(out+"lee_yang_fisher.pdf"); plt.close()
print("ok")
