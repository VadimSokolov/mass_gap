# Figures: langevin_latent.pdf (one plaquette Langevin gaps against the Bakry Emery bound; the law of the latent count M).
import numpy as np, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import betaln
SEED = 2; rng = np.random.default_rng(SEED)   # no random numbers are drawn in this script (latent_mixture.py is exec-ed only up to its first rng line); seed fixed for safety
sys.argv = sys.argv[:1]; sys.path.insert(0, "code")
plt.rcParams.update({"font.family": "serif", "font.size": 8})
g = np.load("code/plaquette_gaps.npy"); fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
ax[0].plot(g[:, 0], g[:, 1], "o-", color="#2166ac", ms=3, label="class functions"); ax[0].plot(g[:, 0], g[:, 2], "s-", color="#1b7837", ms=3, label="sector $l=1$, the gap")
ax[0].plot(g[:, 0], 1-g[:, 0]/2, "--", color="#b2182b", label=r"Bakry \'Emery bound $1-\beta/2$"); ax[0].axhline(0, color="k", lw=0.4); ax[0].set_ylim(-3, 11.5); ax[0].set_xlabel(r"$\beta$"); ax[0].set_ylabel("rate"); ax[0].legend(fontsize=6, frameon=False)
import importlib.util, io, contextlib
spec = importlib.util.spec_from_file_location("lm", "code/latent_mixture.py")
src = open("code/latent_mixture.py").read().split("rng = np.random")[0]; ns = {}; exec(src, ns)
for i, t in enumerate((2.0, 1.0, 0.5, 0.25)):
    p = 1/ns["roots"](t, int(60/t)+10); pm = np.array([1.0])
    for q in p: pm = np.convolve(pm, [1-q, q])
    m = np.arange(pm.size); lw = np.log(np.clip(pm, 1e-300, None))+m*np.log(2)+betaln(m+1.5, 1.5); w = np.exp(lw-lw.max()); w /= w.sum()
    ax[1].plot(m[:30], w[:30], "o-", ms=2.5, lw=0.8, color=plt.cm.viridis(i/3.4), label=r"$\beta=%g$" % (2/t))
ax[1].set_xlabel("latent count $M$"); ax[1].set_ylabel("probability"); ax[1].legend(fontsize=6, frameon=False); fig.tight_layout(); fig.savefig("figures/langevin_latent.pdf")
