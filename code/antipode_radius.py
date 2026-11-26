# Where does the bond moved heat kernel K_t^k stop satisfying the total positivity condition?
# d/dalpha [Delta(f^k)/f^k] = k d/dalpha[Delta f/f] + 2k(k-1) v v',  v = d/dalpha log f, f = K_t from the image sum.
# alpha*(t,k) is the first zero; eta = pi - alpha*. Small t prediction: eta ~ sqrt((k-1) t/(2k)).
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
al = np.linspace(1e-3, np.pi-1e-4, 200001)
def pieces(t):
    kk = np.arange(-5, 6)[:, None]; x = al[None, :]+2*np.pi*kk; e = np.exp(-(x**2-al[None, :]**2)/t)
    g = (x*e).sum(0); g1 = ((1-2*x**2/t)*e).sum(0); g2 = ((-6*x/t+4*x**3/t**2)*e).sum(0)
    gt = (x**3/t**2*e).sum(0); gt1 = ((3*x**2/t**2-2*x**4/t**3)*e).sum(0)
    v = g1/g-1/np.tan(al); vp = g2/g-(g1/g)**2+1/np.sin(al)**2; dq = 4*(gt1*g-gt*g1)/g**2      # d/dalpha (Delta K/K) = 4 d_t d_alpha log K
    return v, vp, dq, g/np.sin(al)
def eta(t, k):
    v, vp, dq, K = pieces(t); c = k*dq+2*k*(k-1)*v*vp; bad = np.flatnonzero(c < 0)
    if bad.size == 0: return 0.0, 0.0
    i = bad[0]; lg = k*(np.log(K)-al**2/t)+2*np.log(np.sin(al)); dens = np.exp(lg-lg.max()); return np.pi-al[i], np.trapezoid(dens[i:], al[i:])/np.trapezoid(dens, al)
print("  t     k=2: eta, prediction, mass beyond      k=4: eta, prediction, mass beyond")
ts = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2]; rows = []
for t in ts:
    r = [t]
    for k in (2, 4): e, m = eta(t, k); r += [e, np.sqrt((k-1)*t/(2*k)), m]
    rows.append(r); print("%5.2f   %.4f %.4f %.2e      %.4f %.4f %.2e" % tuple(r))
assert all((pieces(t)[2] > 0).all() for t in ts)        # the heat kernel itself satisfies the condition everywhere
tt = np.geomspace(0.03, 6, 60); fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
for k, c in [(2, "#16a085"), (4, "#c0392b")]:
    E = np.array([eta(t, k) for t in tt]); ax[0].loglog(tt, E[:, 0], color=c, label=f"$k={k}$"); ax[0].loglog(tt, np.sqrt((k-1)*tt/(2*k)), color=c, ls=":")
    ax[1].loglog(tt, np.clip(E[:, 1], 1e-300, None), color=c, label=f"$k={k}$")
ax[0].set_xlabel("$t$"); ax[0].set_ylabel(r"radius $\eta$ of the bad cap at the antipode"); ax[0].legend(frameon=False, fontsize=7)
ax[1].set_xlabel("$t$"); ax[1].set_ylabel("mass of the bond moved law in the cap"); ax[1].set_ylim(1e-60, 1); ax[1].legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig("figures/antipode_radius.pdf")
