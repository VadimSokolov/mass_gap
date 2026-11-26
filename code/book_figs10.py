# Figure counts.pdf: the exact ladder on counts.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 8})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7)); o = np.load("code/count_orbit_4.npy"); j = np.arange(o.shape[1])
for r in range(0, 8): ax[0].plot(j, o[r], "o-", ms=2.5, lw=0.8, color=plt.cm.viridis(r/7.5), label="rung %d" % r)
ax[0].set_xlabel("count $j$"); ax[0].set_ylabel("probability"); ax[0].legend(fontsize=5.5, frameon=False, ncol=2)
for b, c in ((2, "#b2182b"), (3, "#1b7837"), (4, "#2166ac")):
    o = np.load("code/count_orbit_%d.npy" % b); m = (o*j).sum(1); ax[1].plot(np.arange(m.size), m, "o-", ms=3, color=c, label=r"start $\beta=%d$" % b)
r = np.arange(8); ax[1].plot(r, 6.9-0.5*r, "k:", lw=0.8, label="slope $-1/2$"); ax[1].set_xlabel("rung"); ax[1].set_ylabel("mean count"); ax[1].set_ylim(0, 7.5); ax[1].legend(fontsize=6, frameon=False); fig.tight_layout(); fig.savefig("figures/counts.pdf")
