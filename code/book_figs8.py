# Figures for the chapters on sufficiency, GBC and hitting times: suff_gbc.pdf, hitting_gap.pdf, gbc_gap.pdf.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 8, "font.weight": "normal", "axes.titleweight": "normal"})
d = np.load("code/gbc_exact.npz"); fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
eff = np.array([[1.00, 0.96, 0.53], [0.96, 0.94, 0.57], [0.97, 0.87, 0.56]]); w = 0.26
for j, (lab, c) in enumerate((("sufficient (MLE)", "#1b7837"), ("pseudolikelihood", "#2166ac"), ("one bit summary", "#b2182b"))): ax[0].bar(np.arange(3)+(j-1)*w, eff[:, j], w, color=c, label=lab)
ax[0].axhline(1, color="k", lw=0.6, ls=":"); ax[0].set_xticks(range(3)); ax[0].set_xticklabels([r"$\beta=1$", r"$\beta=2$", r"$\beta=4$"]); ax[0].set_ylabel("efficiency against Cramér Rao"); ax[0].set_ylim(0, 1.6); ax[0].legend(fontsize=6, frameon=False, loc="upper center", ncol=1)
for a, key, c, ttl in ((ax[1], "suff", "#1b7837", "generator on the sufficient statistic"), (ax[2], "lossy", "#b2182b", "generator on the one bit summary")):
    a.plot([0.5, 6], [0.5, 6], color="k", lw=0.6)
    for q, al in ((0, 0.35), (2, 0.9), (4, 0.35)): a.scatter(d["ex"][:, q], d[key][:, q], s=3, color=c, alpha=al, lw=0)
    a.set_xlabel("exact posterior quantile"); a.set_ylabel("generated quantile"); a.set_title(ttl, fontsize=7)
fig.tight_layout(); fig.savefig("figures/suff_gbc.pdf"); plt.close(fig)
h = np.load("code/hitting_gap.npz"); t = h["t"]; N = int(h["N"]); fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
ex = sum((16/np.pi)*(-1)**(k+1)*k*k/(4*k*k-1)*np.exp(-(4*k*k-1)*t/2) for k in range(1, 40)); ok = h["surv"] > 0
ax[0].semilogy(t[ok], h["surv"][ok], color="#2166ac", lw=1.4, label="simulated survival of the hitting time"); ax[0].semilogy(t[t > 0.15], ex[t > 0.15], color="#b2182b", lw=0.9, ls="--", label="exact Dirichlet series")
ax[0].semilogy(t, 16/(3*np.pi)*np.exp(-1.5*t), color="k", lw=0.5, ls=":", label=r"$\frac{16}{3\pi}e^{-3t/2}$"); ax[0].set_xlabel("time $t$"); ax[0].set_ylabel(r"$P(T>t)$"); ax[0].set_ylim(1e-3, 1.5); ax[0].legend(fontsize=6, frameon=False)
c = h["corr"]; se = 0.5/np.sqrt(N); ax[1].fill_between(t, np.clip(c-2*se, 1e-5, None), c+2*se, color="#f4a582", alpha=0.6, lw=0); ax[1].semilogy(t, np.clip(c, 1e-5, None), color="#b2182b", lw=1.0, label=r"correlator $E\,a_0(t)$ with two standard errors")
s = h["surv"]; ses = np.sqrt(s*(1-s)/N); ax[1].fill_between(t, np.clip(s-2*ses, 1e-5, None), s+2*ses, color="#92c5de", alpha=0.7, lw=0); ax[1].semilogy(t, np.clip(s, 1e-5, None), color="#2166ac", lw=1.0, label="survival with two standard errors")
ax[1].set_ylim(1e-3, 2); ax[1].set_xlabel("time $t$"); ax[1].legend(fontsize=6, frameon=False); fig.tight_layout(); fig.savefig("figures/hitting_gap.pdf"); plt.close(fig)
print("survival at t=2: sim %.4f exact %.4f" % (h["surv"][np.argmin(abs(t-2))], ex[np.argmin(abs(t-2))]))
try:
    g = np.load("code/gbc_gap.npz"); i = np.argsort(g["mt"])[::8]; fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
    ax[0].vlines(g["mt"][i], g["lo"][i], g["hi"][i], color="#92c5de", lw=0.5); ax[0].scatter(g["mt"][i], g["md"][i], s=4, color="#2166ac", zorder=3); ax[0].plot([0.2, 1.2], [0.2, 1.2], "k", lw=0.6); ax[0].set_xlabel("true gap $m$"); ax[0].set_ylabel("GBC median and 90% interval")
    ax[1].scatter(g["mt"][i], g["e3"][i], s=4, color="#b2182b"); ax[1].plot([0.2, 1.2], [0.2, 1.2], "k", lw=0.6); ax[1].set_ylim(0, 2.2); ax[1].set_xlabel("true gap $m$"); ax[1].set_ylabel(r"effective mass at $\tau=3$"); fig.tight_layout(); fig.savefig("figures/gbc_gap.pdf")
except FileNotFoundError: print("gbc_gap not ready")
