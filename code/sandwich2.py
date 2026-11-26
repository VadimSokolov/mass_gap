# The sandwich in three and four dimensions and over two steps: plaquette -> 2x2 -> 4x4 loop laws.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
bins = np.linspace(0, np.pi, 61); th = 0.5*(bins[1:]+bins[:-1]); dth = bins[1]-bins[0]; n = np.arange(1, 40)
CHI = np.sin(np.outer(n, th))/np.sin(th); haar = (2/np.pi)*np.sin(th)**2*dth
coeff = lambda h: (lambda c: c/c[0])((CHI@(h/h.sum()))/n)
def mk(u, k):
    f = np.clip((n*u)@CHI, 1e-12, None); c = (CHI@(f**k*haar))/n; return (c/c[0])**4
dens = lambda u: np.clip((n*u)@CHI, 1e-12, None)*(2/np.pi)*np.sin(th)**2
def share(num, den, wts):
    r = num/np.where(den > 1e-9, den, np.nan); ok = np.isfinite(r[1:]) & np.isfinite(r[:-1]); w = np.minimum(wts[1:], wts[:-1])[ok]
    return (w*(np.diff(r)[ok] > 0)).sum()/w.sum()
rows = []
for D, betas in [(4, ["2.5", "3.0", "4.0", "6.0"]), (3, ["4.0", "6.0", "9.0", "14.0"])]:
    k = 2**(D-2)
    for b in betas:
        h1, h2, h4 = np.load(f"/tmp/loopsd_{D}_{b}.npy"); u1, u2, u4 = coeff(h1), coeff(h2), coeff(h4)
        for name, src, tgt, h, steps in [("1 to 2", u1, u2, h2, 1), ("2 to 4", u2, u4, h4, 1), ("1 to 4", u1, u4, h4, 2)]:
            ind = src**(4**steps); m = mk(src, k) if steps == 1 else mk(mk(src, k), k)
            lam = (np.log(tgt[1])-np.log(ind[1]))/(np.log(m[1])-np.log(ind[1])); p = h/h.sum()/dth
            rows.append((D, float(b), name, ind[1], tgt[1], m[1], lam, share(dens(ind), p, p), share(p, dens(m), p)))
            print("D=%d beta %5s %s  ind %.4f meas %.4f MK %.4f  position %.3f  lr shares %.2f %.2f" % ((D, b, name)+rows[-1][3:]))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True); cols = {"1 to 2": "#1f4e79", "2 to 4": "#c0392b", "1 to 4": "#16a085"}
for a, D in zip(ax, (4, 3)):
    for name in cols:
        r = [x for x in rows if x[0] == D and x[2] == name]; a.plot([x[1] for x in r], [x[6] for x in r], "o-", color=cols[name], label=name)
    a.axhline(0, color="k", lw=0.6); a.axhline(1, color="k", lw=0.6); a.set_ylim(-0.1, 1.1); a.set_xlabel(r"$\beta$"); a.set_title(f"{D} dimensions", fontsize=9)
    a.text(a.get_xlim()[0], 0.02, " independent", fontsize=7); a.text(a.get_xlim()[0], 0.93, " Migdal Kadanoff", fontsize=7)
ax[0].set_ylabel("position of the truth in the bracket"); ax[0].legend(frameon=False, fontsize=7, loc="center right")
fig.tight_layout(); fig.savefig("figures/sandwich2.pdf")
