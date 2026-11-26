# Local heat kernel time of a class function law f on SU(2): t_f(theta) solves score_{K_t}(theta) = score_f(theta).
# It is well defined because the score of the heat kernel is strictly increasing in t (the product formula theorem).
# K_{t-} <=_lr f <=_lr K_{t+} with t- = min t_f, t+ = max t_f. We follow the bracket along the ladder from a Wilson start.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ive
M, NMAX = 4000, 400; x, w = np.polynomial.legendre.leggauss(M); th = (x+1)*np.pi/2; w = w*np.pi/2; haar = (2/np.pi)*np.sin(th)**2*w
n = np.arange(1, NMAX+1); CHI = np.sin(np.outer(n, th))/np.sin(th); sel = (th > 0.08) & (th < np.pi-0.08); ths = th[sel]
DCHI = (n[:, None]*np.cos(np.outer(n, th))*np.sin(th)-np.sin(np.outer(n, th))*np.cos(th))/np.sin(th)**2
def score(u): return (((n*u)@DCHI)/((n*u)@CHI))[sel]          # d/dtheta log f for f = sum n u_n chi_n
def score_heat(t):                                            # from the image sum, stable at small t
    k = np.arange(-4, 5)[:, None]; xk = ths[None, :]+2*np.pi*k; e = np.exp(-(xk**2-ths[None, :]**2)/t)
    return ((1-2*xk**2/t)*e).sum(0)/(xk*e).sum(0)-1/np.tan(ths)
tg = np.geomspace(0.02, 60, 500); S = np.array([score_heat(t) for t in tg]); assert (np.diff(S, axis=0) > 0).all()
def local_time(u):
    s = score(u); return np.array([np.interp(s[i], S[:, i], tg) for i in range(ths.size)])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1)); cols = plt.cm.viridis(np.linspace(0.05, 0.9, 7))
for beta, ls in [(2.5, "-")]:
    u = np.clip(ive(n, beta)/ive(1, beta), 1e-300, None); out = []
    for m in range(7):
        tl = local_time(u); tm = -4*np.log(u[1])/3; out.append((m, tl.min(), tm, tl.max()))
        ax[0].plot(ths, tl/tm, color=cols[m], label=f"rung {m}")
        f = np.clip((n*u)@CHI, 1e-300, None); c = (CHI@(f**4*haar))/n; u = np.clip((c/c[0]), 1e-300, None)**4
    out = np.array(out)
    for r in out: print("Wilson beta %.1f rung %d: t- %.4f  t(u2) %.4f  t+ %.4f  width t+/t- %.3f" % (beta, r[0], r[1], r[2], r[3], r[3]/r[1]))
ax[0].axhline(1, color="k", lw=0.6); ax[0].set_xlabel(r"$\theta$"); ax[0].set_ylabel("local heat kernel time, relative"); ax[0].legend(frameon=False, fontsize=6, ncol=2)
ax[1].semilogy(out[:, 0], out[:, 3]/out[:, 1]-1, "o-", color="#c0392b"); ax[1].set_xlabel("rung"); ax[1].set_ylabel(r"bracket width $t_+/t_--1$")
fig.tight_layout(); fig.savefig("figures/local_time.pdf")
