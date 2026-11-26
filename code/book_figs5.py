# Colour figures for the slice, Besag and progress chapters.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
import sys; sys.path.insert(0, "code")
from mk_verify import make_R
F = "figures/"
def logf(t, th, K=6):
    k = np.arange(-K, K+1)[:, None, None]; x = th[None] + 2*np.pi*k
    return np.log(np.abs((x*np.exp(-x**2/t)).sum(0)))
# 1. mixed partial heat map, book time t (image sum exponent x^2/t)
t = np.geomspace(0.04, 8, 220); th = np.linspace(0.01, np.pi-0.01, 300); T, TH = np.meshgrid(t, th, indexing="ij")
h = 1e-4*T; d = 1e-4
m = ((logf(T+h, TH+d)-logf(T-h, TH+d))-(logf(T+h, TH-d)-logf(T-h, TH-d)))/(4*h*d)/np.sin(TH)
fig, ax = plt.subplots(figsize=(5.6, 3.6))
pc = ax.pcolormesh(TH, T, np.log10(np.clip(m, 1e-6, None)), cmap="viridis", shading="auto")
ax.set_yscale("log"); ax.axhline(3.214, color="w", ls="--", lw=1); ax.text(0.15, 3.7, r"theorem: $t\geq3.214$", color="w", fontsize=8)
ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$t$"); fig.colorbar(pc, label=r"$\log_{10}\,\partial_t\partial_\theta\log K_t/\sin\theta$")
fig.tight_layout(); fig.savefig(F+"mlr_heatmap.pdf"); plt.close()
# 2. likelihood ratios K_t/K_1
th = np.linspace(0.001, np.pi-0.001, 400)
def K(t): n = np.arange(1, 200)[:, None]; return (n*np.exp(-(n**2-1)*t/4)*np.sin(n*th)/np.sin(th)).sum(0)
fig, ax = plt.subplots(figsize=(5.4, 3.4)); ts = [1.2, 1.5, 2, 3, 5, 8]; col = cm.plasma(np.linspace(0.1, 0.85, len(ts)))
for tt, c in zip(ts, col): r = K(tt)/K(1.0); ax.semilogy(th, r/r[0], color=c, label=f"$t={tt}$")
ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$K_t(\theta)/K_1(\theta)$, normalised"); ax.legend(frameon=False, fontsize=7, ncol=2)
ax.set_title("Monotone likelihood ratios of the SU(2) heat kernel", fontsize=9); fig.tight_layout(); fig.savefig(F+"mlr_ratios.pdf"); plt.close()
# 3. decrement against perturbation theory
R = make_R(3000, 600); bb = np.geomspace(0.3, 60, 80); dl = np.array([b-R(b) for b in bb])
fig, ax = plt.subplots(figsize=(5.4, 3.4))
ax.semilogx(bb, dl, color="#1f4e79", lw=2, label=r"ladder $\beta-R(\beta)$")
ax.semilogx(bb, 0.25+0.021/bb, color="#2e86c1", ls=":", label=r"$0.25+0.021/\beta$")
ax.semilogx(bb, 0.2575+0.0403/bb, color="#e67e22", label=r"two loop $0.2575+0.0403/\beta$")
ax.axhline(0.2575, color="#7f8c8d", ls="--", lw=1, label="one loop")
ax.set_ylim(0.2, 0.4); ax.set_xlabel(r"$\beta$"); ax.set_ylabel("decrement per doubling"); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(F+"decrement_vs_pert.pdf"); plt.close()
# 4. plaquette: our 4D data, 2D exact, weak coupling
from scipy.special import iv
b4 = np.array([2.2, 2.3, 2.4, 2.5, 2.635]); P4 = np.array([0.569972, 0.603046, 0.629759, 0.652635, 0.676289])
x = np.linspace(0.5, 4, 200)
fig, ax = plt.subplots(figsize=(5.4, 3.4))
ax.plot(x, iv(2, x)/iv(1, x), color="#16a085", label="2D exact, single plaquette")
ax.plot(x, 1-3/(4*x), color="#e67e22", ls="--", label=r"4D weak coupling $1-3/(4\beta)$")
ax.plot(x, x/4, color="#8e44ad", ls=":", label=r"strong coupling $\beta/4$")
ax.plot(b4, P4, "o", color="#c0392b", label=r"4D heat bath, $8^4$")
ax.set_ylim(0, 1); ax.set_xlabel(r"$\beta$"); ax.set_ylabel("mean plaquette"); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(F+"plaquette_4d.pdf"); plt.close()
# 5. pseudolikelihood score curves, synthetic from exact conditionals
from scipy.special import ive
rng = np.random.default_rng(3); fig, ax = plt.subplots(figsize=(5.4, 3.4)); col = cm.viridis(np.linspace(0.1, 0.9, 4))
for b, c in zip([1, 2, 4, 8], col):
    k = rng.uniform(0.6, 2.0, 288); kap = b*k; a0 = np.empty_like(kap)
    for i, kp in enumerate(kap):
        while True:
            y = rng.uniform(-1, 1)
            if rng.random() < np.sqrt(1-y*y)*np.exp(kp*(y-1)): a0[i] = y; break
    g = np.linspace(0.3, 12, 300); sc = [np.sum(k*(a0-ive(2, gg*k)/ive(1, gg*k)))/len(k) for gg in g]
    ax.plot(g, sc, color=c, label=fr"true $\beta={b}$"); ax.axvline(b, color=c, ls=":", lw=1)
ax.axhline(0, color="k", lw=0.6); ax.set_xlabel(r"$\beta$"); ax.set_ylabel("pseudolikelihood score per link"); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(F+"pl_score.pdf"); plt.close()
# 6. slice caps on the sphere (two dimensional picture)
def cap_mask(X, Y, Z, c, a): return (X*c[0]+Y*c[1]+Z*c[2]) > a
u, v = np.meshgrid(np.linspace(0, 2*np.pi, 400), np.linspace(0, np.pi, 200)); X, Y, Z = np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)
cs = [np.array(c)/np.linalg.norm(c) for c in [(0.3, 0.2, 1), (-0.2, 0.5, 1), (0.4, -0.4, 1), (-0.5, -0.3, 1)]]
lv = [0.55, 0.5, 0.6, 0.45]; inter = np.ones_like(X, bool); count = np.zeros_like(X)
for c, a in zip(cs, lv): mk = cap_mask(X, Y, Z, c, a); inter &= mk; count += mk
fig = plt.figure(figsize=(5.2, 3.8)); ax = fig.add_subplot(projection="3d")
colors = cm.Blues(0.15+0.2*count); colors[inter] = (0.85, 0.2, 0.15, 1)
ax.plot_surface(X, Y, Z, facecolors=colors, rstride=2, cstride=2, linewidth=0, shade=False)
ax.view_init(35, 30); ax.set_axis_off(); ax.set_box_aspect((1, 1, 1))
ax.set_title("Slice caps and their convex intersection", fontsize=9)
fig.tight_layout(); fig.savefig(F+"slice_caps.pdf", dpi=150); plt.close()
print("done")
