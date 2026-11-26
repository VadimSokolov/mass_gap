# Figures for the chapters on counts, Kingman's coalescent and real zeros: counts_urn.pdf, kingman_fig.pdf, zeros_fig.pdf.
import numpy as np, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from math import comb
from scipy.special import jn_zeros
from scipy.optimize import brentq
sys.path.insert(0, "code"); import kingman_dual as KD
plt.rcParams.update({"font.family": "serif", "font.size": 8}); cm = plt.cm.viridis
# --- counts_urn
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
for i, (a, b) in enumerate(((4, 16), (8, 8), (8, 16), (16, 16), (16, 40))):
    q = [comb(a, j)*comb(b+2, b-j)/comb(a+b+2, b) for j in range(min(a, b)+1)]; ax[0].plot(range(len(q)), q, "o-", ms=2.5, lw=0.8, color=cm(i/4.5), label="$a=%d,\\ b=%d$" % (a, b)); ax[0].axvline(a*b/(a+b+2), color=cm(i/4.5), lw=0.4, ls=":")
ax[0].set_xlabel("size $j$ of the intersection"); ax[0].set_ylabel("$Q(j\\mid a,b)$"); ax[0].legend(fontsize=6, frameon=False)
r = np.linspace(0, 1.15, 50); ax[1].plot(r, -13/16+15*r/16, color="k", lw=0.9, label=r"$-\frac{13}{16}+\frac{15}{16}r$")
for (rr, d, lab, c, mk) in ((0, -0.7976, "point masses", "#2166ac", "s"), (0, -0.8075, None, "#2166ac", "s"), (0.388, -0.5074, "heat kernels", "#1b7837", "o"), (0.360, -0.5042, None, "#1b7837", "o"), (0.346, -0.5022, None, "#1b7837", "o"), (1.095, 0.0912, "Wilson laws", "#b2182b", "^"), (1.024, 0.1178, None, "#b2182b", "^")):
    ax[1].plot(rr, d, mk, color=c, ms=5, label=lab)
ax[1].axhline(0, color="gray", lw=0.4); ax[1].axvline(1/3, color="gray", lw=0.4, ls=":"); ax[1].set_xlabel("$r=$ variance / mean of the count"); ax[1].set_ylabel("change of the mean count in one rung"); ax[1].legend(fontsize=6, frameon=False); fig.tight_layout(); fig.savefig("figures/counts_urn.pdf"); plt.close(fig)
# --- kingman_fig
rng = np.random.default_rng(6); fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7)); N0 = 400
for p in range(6):
    m = N0; tau = 0.0; ts, ms = [], []
    while m > 0: tau += rng.exponential(2/(m*(m+2))); m -= 1; ts.append(2*tau); ms.append(m)
    ax[0].step(ts, ms, where="post", lw=0.7, color=cm(p/5.5))
t = np.geomspace(0.02, 12, 200); ax[0].plot(t, 4/t-7/6, "k--", lw=0.8, label=r"$4/t-\frac{7}{6}$"); ax[0].set_xscale("log"); ax[0].set_yscale("symlog", linthresh=1); ax[0].set_xlim(0.02, 12); ax[0].set_ylim(-0.2, 300); ax[0].set_xlabel("heat kernel time $t$"); ax[0].set_ylabel("number of lineages"); ax[0].legend(fontsize=6, frameon=False)
tt = np.linspace(0.6, 12, 60); p0 = np.array([KD.tavare(x)[0] for x in tt]); ax[1].semilogy(tt, 1-p0, color="#2166ac", lw=1.3, label="$1-$ weight of Haar measure in $K_t$"); ax[1].semilogy(tt, 4*np.exp(-3*tt/4), "k:", lw=0.8, label=r"$4\,e^{-3t/4}$"); ax[1].set_xlabel("$t$"); ax[1].legend(fontsize=6, frameon=False); fig.tight_layout(); fig.savefig("figures/kingman_fig.pdf"); plt.close(fig)
print("1 - P(A = 0) over 3 exp(-3t/4) at t = 8, 12:", (1-p0[np.argmin(abs(tt-8))])/(3*np.exp(-6)), (1-p0[-1])/(3*np.exp(-9)))
# --- zeros_fig
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7)); ts = np.geomspace(0.15, 8, 40)
G = []
for t in ts:
    N = int(60/t)+60; C = np.cosh((2*np.arange(1, N+1)-1)*t/4); f = lambda g: np.sum(1/(C-g)); G.append([brentq(f, C[m]*(1+1e-13)+1e-13, C[m+1]*(1-1e-13)-1e-13) for m in range(6)])
G = np.array(G)
for m in range(6): ax[0].loglog(ts, G[:, m]-1, color=cm(m/5.5), lw=1.1, label="$m=%d$" % (m+1))
ax[0].set_xlabel("heat kernel time $t$"); ax[0].set_ylabel(r"$\gamma_m(t)-1$"); ax[0].legend(fontsize=6, frameon=False, ncol=2)
kap = np.geomspace(0.5, 20, 60); jz = jn_zeros(1, 6)
for m in range(6): ax[1].loglog(kap, jz[m]**2/(2*kap**2), color=cm(m/5.5), lw=1.1)
ax[1].set_xlabel(r"$\kappa$ (two Wilson laws, $\kappa=\lambda$)"); ax[1].set_ylabel(r"$\gamma_m-1=j_{1,m}^2/2\kappa^2$"); fig.tight_layout(); fig.savefig("figures/zeros_fig.pdf")
