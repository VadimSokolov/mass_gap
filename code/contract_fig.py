# Figure for the contraction theorem: zeta(t, x) = -d/dt[t d/dx log K_t] from the product formula, and the likelihood ratio K_{s/k}/K_s^k.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 8})
def parts(t, x, N=None):
    N = N or int(60/t)+5; n = np.arange(1, N+1)[:, None]; y = (2*n-1)*t/4; C = np.cosh(y); Cp = (2*n-1)/4*np.sinh(y); a = 1/(C+x)
    S1, S2, D2, D3 = a.sum(0), (a*a).sum(0), (Cp*a*a).sum(0), (Cp*a**3).sum(0); return S1-S2/S1, t*(D2-2*D3/S1+D2*S2/S1**2)-S1+S2/S1
th = np.linspace(0.01, np.pi-0.01, 600); x = np.cos(th); fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7)); cm = plt.cm.viridis
for i, t in enumerate((0.05, 0.3, 1, 2, 4, 8)):
    G, ze = parts(t, x); ax[0].semilogy(th, ze*np.sin(th), color=cm(i/5.5), lw=1.1, label="$t=%g$" % t)
ax[0].semilogy(th, 1/th-1/np.tan(th), "k:", lw=0.8, label=r"$1/\theta-\cot\theta$"); ax[0].set_xlabel(r"$\theta$"); ax[0].set_ylabel(r"$z(t,\theta)$"); ax[0].set_ylim(1e-4, 1e2); ax[0].legend(fontsize=6, frameon=False, ncol=2)
for i, s in enumerate((0.5, 1, 2, 4, 8)):
    lr = np.cumsum((parts(s/4, x)[0]-4*parts(s, x)[0])[::-1]*np.gradient(x[::-1]))[::-1]; ax[1].plot(x, lr-lr[-1], color=cm(i/4.5), lw=1.1, label=r"$\beta=%g$" % (2/s))
ax[1].set_xlabel(r"$\cos\theta$"); ax[1].set_ylabel(r"$\log\,K_{s/4}/K_s^4$, shifted"); ax[1].legend(fontsize=6, frameon=False); fig.tight_layout(); fig.savefig("figures/contract.pdf")
