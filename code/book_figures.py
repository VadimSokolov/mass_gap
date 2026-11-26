"""Additional figures for the book. Polson and Sokolov, Nov 2026."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ive
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
out = "../figures/"
th = np.linspace(1e-4, np.pi-1e-4, 800); n = np.arange(1, 600)[:, None]
def K(t): return ((n*np.exp(-t*(n**2-1)/4))*np.sin(n*th[None, :])).sum(0)/np.sin(th)
haar = (2/np.pi)*np.sin(th)**2
# 1. heat kernel laws of the rotation angle
fig, ax = plt.subplots(figsize=(5.6, 3.2)); cm = plt.cm.viridis
for i, t in enumerate([0.05, 0.15, 0.4, 1.0, 2.5, 6.0]):
    ax.plot(th, K(t)*haar, color=cm(i/5), label="$t=%g$" % t)
ax.plot(th, haar, "k:", label="Haar")
ax.set_xlabel("rotation angle $\\theta$"); ax.set_ylabel("density"); ax.legend(frameon=False, fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(out+"heat_kernel_laws.pdf"); plt.close()
# 2. curvature factor and the BPY Laplace transform
x = np.linspace(0, 3, 300)
fig, ax = plt.subplots(1, 2, figsize=(6.2, 2.9))
tt = np.linspace(0, 2.8, 300)
ax[0].plot(tt, tt/np.sin(np.maximum(tt, 1e-9)), "C0", label="$\\theta/\\sin\\theta$")
ax[0].plot(tt, np.exp(tt**2/6), "C3--", label="$e^{\\theta^2/6}$")
ax[0].set_ylim(0.9, 4); ax[0].set_xlabel("$\\theta$"); ax[0].legend(frameon=False, fontsize=8)
lam = np.linspace(0, 6, 300); r = np.sqrt(2*lam)
ax[1].plot(lam, np.where(r > 0, r/np.sinh(np.maximum(r, 1e-9)), 1), "C0", label="$\\sqrt{2\\lambda}/\\sinh\\sqrt{2\\lambda}$")
ax[1].plot(lam, np.exp(-lam/3), "C3--", label="$e^{-\\lambda/3}$")
ax[1].set_xlabel("$\\lambda$"); ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"curvature_bpy.pdf"); plt.close()
# 3. density of S by inversion of the series: f(x) = sum_n (-1)^{n+1} ... use Thorin atom series via numerical Laplace inversion (Talbot)
import mpmath as mp
phi = lambda s: mp.sqrt(2*s)/mp.sinh(mp.sqrt(2*s))
xs = np.linspace(0.02, 1.6, 120)
dens = [float(mp.invertlaplace(phi, float(v), method="talbot")) for v in xs]
fig, ax = plt.subplots(figsize=(5.2, 2.8))
ax.plot(xs, dens, "C0"); ax.axvline(1/3, color="k", ls=":", label="$\\mathbb{E}S=1/3$")
ax.set_xlabel("$x$"); ax.set_ylabel("density of $S$"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"bpy_density.pdf"); plt.close()
print("mass of S density:", np.trapezoid(dens, xs))
# 4. Thorin levels u_n along an MK trajectory
import sys; sys.path.insert(0, ".")

M = 4000; xg, wg = np.polynomial.legendre.leggauss(M); T = (xg+1)*np.pi/2; wg = wg*np.pi/2
H = (2/np.pi)*np.sin(T)**2*wg; NN = np.arange(1, 400); CH = np.sin(np.outer(NN, T))/np.sin(T)
u = ive(NN, 6.0)/ive(1, 6.0); traj = [u.copy()]
for s in range(26):
    f = np.clip((NN*u)@CH, 1e-300, None); lf = 4*np.log(f); g = np.exp(lf-lf.max())
    c = (CH@(g*H))/NN; c = np.clip(c/c[0], 0, None); u = c**4; traj.append(u.copy())
traj = np.array(traj)
fig, ax = plt.subplots(figsize=(5.6, 3.2))
for j, nn in enumerate([2, 3, 4, 5, 6]):
    ax.semilogy(np.arange(len(traj)), np.maximum(traj[:, nn-1], 1e-300), color=plt.cm.plasma(j/5), label="$u_%d$" % nn)
ax.set_ylim(1e-12, 1.5); ax.set_xlabel("MK step (lattice spacing doubles)"); ax.set_ylabel("Thorin level $u_n$")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"mk_levels.pdf"); plt.close()
# 5. SU(3) decrement convergence (values from su3_mk.py)
b = np.array([4, 8, 12, 16, 20, 24]); d = np.array([0.9605, 0.6343, 0.5808, 0.5718, 0.5689, 0.5676])
fig, ax = plt.subplots(figsize=(5.2, 2.9))
ax.plot(b, d, "C0o-", label="SU(3) MK decrement")
ax.axhline(9/16, color="k", ls=":", label="$(1-b^{-2})N^2\\mathbb{E}S/4=9/16$")
ax.axhline(11*9*np.log(2)/(12*np.pi**2), color="C3", ls="--", label="one loop $0.5794$")
ax.set_xlabel("$\\beta$"); ax.set_ylim(0.5, 1.0); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(out+"su3_decrement.pdf"); plt.close()
print("ok")
