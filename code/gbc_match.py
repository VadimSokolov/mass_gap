# GBC for renormalization: learn the map from loop summaries of a 4^4 lattice to the Wilson coupling that produced them,
# as a quantile generator trained on simulated pairs, then apply it to blocked 8^4 ensembles. This is two lattice
# matching in the sense of Monte Carlo renormalization, amortised.
import numpy as np, glob, os, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingRegressor
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 0   # default seed 0 reproduces the original run; data in $GBC_DIR (default /tmp)
DIR = os.environ.get("GBC_DIR", "/tmp")
T = np.concatenate([np.load(f) for f in sorted(glob.glob(DIR+"/gbc_train_*.npy"))]); b, X = T[:, 0], T[:, 1:]
rng = np.random.default_rng(0+SEED); idx = rng.permutation(len(b)); tr, ca, te = idx[:300], idx[300:360], idx[360:]
feat = lambda X: np.column_stack([X, -np.log(np.clip(X, 1e-4, None))])
G = GradientBoostingRegressor(loss="quantile", alpha=0.5, n_estimators=400, max_depth=3, learning_rate=0.05, subsample=0.8, random_state=1+SEED).fit(feat(X[tr]), b[tr])
# split conformal calibration of the interval, with scale proportional to the square of the prediction
pc = G.predict(feat(X[ca])); score = np.sort(np.abs(b[ca]-pc)/pc**2); qhat = score[int(np.ceil(0.9*(len(ca)+1)))-1]
Q = {0.5: G}; lo = lambda x: G.predict(x)*(1-qhat*G.predict(x)); hi = lambda x: G.predict(x)*(1+qhat*G.predict(x))
pr = {0.5: G.predict(feat(X[te])), 0.05: lo(feat(X[te])), 0.95: hi(feat(X[te]))}
print("held out: n %d  median abs error %.3f  coverage of the conformal 90%% interval %.2f  qhat %.4f" % (len(te), np.median(np.abs(pr[0.5]-b[te])), np.mean((b[te] >= pr[0.05]) & (b[te] <= pr[0.95])), qhat))
# single observable matching: invert the smoothed training curve of each loop
o = np.argsort(b); single = []
for j in range(3):
    co = np.polyfit(b[o], -np.log(np.clip(X[o, j], 1e-4, None)), 9); grid = np.linspace(0.3, 9, 4000); single.append((grid, np.polyval(co, grid)))
inv = lambda j, w: single[j][0][np.argmin(np.abs(single[j][1]+np.log(max(w, 1e-4))))]
rows = []
for f in sorted(glob.glob(DIR+"/gbc_block_*.npy")):
    for be, s, w11, w12, w22 in np.load(f):
        x = feat(np.array([[w11, w12, w22]])); m = [inv(0, w11), inv(1, w12), inv(2, w22)]
        rows.append((be, s, G.predict(x)[0], lo(x)[0], hi(x)[0], m[0], m[1], m[2], max(m)-min(m)))
rows = np.array(rows); np.save(DIR+"/gbc_rows.npy", rows)
print(" beta smear | GBC median [5%, 95%] | matched by W11, W12, W22 | spread")
for r in rows: print(" %4.1f %4.1f | %6.3f [%6.3f, %6.3f] | %6.3f %6.3f %6.3f | %.3f" % tuple(r))
if SEED != 0: sys.exit()   # figures only for the original seed
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
sc = ax[0].scatter(b[te], pr[0.5], c=X[te, 0], cmap="viridis", s=14); ax[0].plot([0, 9], [0, 9], "k--", lw=0.8)
ax[0].vlines(b[te], pr[0.05], pr[0.95], color="#7f8c8d", lw=0.5, zorder=0); ax[0].set_xlabel(r"true $\beta$"); ax[0].set_ylabel("GBC posterior median and 90% interval"); fig.colorbar(sc, ax=ax[0], label="$W_{11}$")
for be, c in zip(sorted(set(rows[:, 0])), plt.cm.plasma(np.linspace(0.05, 0.85, 5))):
    r = rows[rows[:, 0] == be]; ax[1].plot(r[:, 1], r[:, 5]-be, "o-", color=c, ms=3, label=fr"$\beta={be:g}$"); ax[1].plot(r[:, 1], r[:, 7]-be, "s--", color=c, ms=3)
ax[1].axhline(-0.26, color="k", lw=0.8); ax[1].set_ylim(-4, 6); ax[1].set_xlabel("smearing weight"); ax[1].set_ylabel(r"matched $\beta'-\beta$"); ax[1].legend(frameon=False, fontsize=6)
fig.tight_layout(); fig.savefig("figures/gbc_match.pdf")
