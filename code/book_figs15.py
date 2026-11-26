# Figure: where the lattice laws sit between the heat kernel (0) and the Wilson law (1) with the same fundamental coefficient.
# Data: output of count_lattice.py (8^4, forty configurations) and count_blocked.py.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 9, "mathtext.fontset": "cm"})
b = np.array([1.5, 2.0, 2.3, 2.5, 2.7, 3.0, 4.0])
p = np.array([[2.1727, 2.1720, 1.9013, .0039], [3.6662, 3.6862, 3.2524, .0064], [5.2341, 5.3484, 4.8091, .0109], [6.2989, 6.4766, 5.8963, .0114], [7.1925, 7.4033, 6.8009, .0067], [8.4964, 8.7307, 8.1082, .0125], [12.6360, 12.9152, 12.2674, .0112]])
l = np.array([[.0714, .0766, .0753, .0016], [.3023, .3159, .2989, .0028], [.7938, .8656, .7809, .0063], [1.2268, 1.3643, 1.2068, .0059], [1.5569, 1.7557, 1.5417, .0042], [2.0380, 2.3116, 2.0229, .0071], [3.5668, 4.0043, 3.5450, .0084]])
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7))
for d, mk, lab in ((p, "ko-", "plaquette"), (l[1:], "ks--", r"$2\times2$ loop")):
    bb = b if d is p else b[1:]; ax[0].errorbar(bb, (d[:, 0]-d[:, 2])/(d[:, 1]-d[:, 2]), d[:, 3]/(d[:, 1]-d[:, 2]), fmt=mk, ms=3, lw=0.8, mfc="w", label=lab)
ax[0].axhline(1, color="0.6", lw=0.6); ax[0].axhline(0, color="0.6", lw=0.6); ax[0].text(3.3, 1.03, "Wilson", fontsize=8); ax[0].text(3.1, 0.2, "heat kernel", fontsize=8); ax[0].set_ylim(-0.25, 1.2)
ax[0].set_xlabel(r"$\beta$"); ax[0].set_ylabel("position between the two shapes"); ax[0].legend(frameon=False, fontsize=8, loc="lower left", ncol=2)
bb = np.array([2.0, 2.3, 2.7]); fine = [3.64, 5.17, 7.25]; lad = [3.61, 5.33, 7.47]; s0 = [0.27, 0.75, 1.48]; s5 = [1.50, 3.77, 8.22]; s2 = [1.29, 3.47, 8.34]
ax[1].plot(bb, fine, "k^-", ms=4, lw=0.8, mfc="w", label="fine"); ax[1].plot(bb, lad, "k-", lw=1.4, label="ladder, one rung"); ax[1].plot(bb, s5, "ko--", ms=3, lw=0.8, label="blocked, smear 0.5"); ax[1].plot(bb, s2, "ks:", ms=3, lw=0.8, label="blocked, smear 2"); ax[1].plot(bb, s0, "kx-.", ms=4, lw=0.8, label="blocked, plain")
ax[1].set_xlabel(r"$\beta$"); ax[1].set_ylabel("count"); ax[1].legend(frameon=False, fontsize=7)
plt.tight_layout(); plt.savefig("../figures/count_lattice.pdf")
