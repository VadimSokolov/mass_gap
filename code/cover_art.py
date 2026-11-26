"""Cover art: laws of Brownian motion on SU(2) (heat kernel times Haar) at many times, glowing on navy."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
th = np.linspace(1e-4, np.pi-1e-4, 1200); n = np.arange(1, 900)[:, None]
haar = (2/np.pi)*np.sin(th)**2
def K(t): return ((n*np.exp(-t*(n**2-1)/4))*np.sin(n*th[None, :])).sum(0)/np.sin(th)
ts = np.geomspace(0.03, 6, 60)
fig = plt.figure(figsize=(6, 4.2), dpi=400); ax = fig.add_axes([0, 0, 1, 1])
fig.patch.set_alpha(0); ax.set_facecolor((0, 0, 0, 0))
cm = plt.cm.plasma
for i, t in enumerate(ts):
    y = K(t)*haar; c = cm(0.15+0.8*i/len(ts))
    for lw, al in [(6, 0.04), (3, 0.08), (1.0, 0.9)]:
        ax.plot(th, y + 0.0*i, color=c, lw=lw, alpha=al, solid_capstyle="round")
ax.plot(th, haar, color=(1, 0.85, 0.5), lw=1.2, ls=(0, (2, 3)), alpha=0.9)
ax.set_xlim(0, np.pi); ax.set_ylim(-0.05, 3.3); ax.axis("off")
fig.savefig("../cover/cover_art.png", transparent=True)
print("ok")
