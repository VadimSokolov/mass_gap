# Figure for the Schur and Szego criterion: an antiderivative of f = (1+x)(1+x/10)(1+x/10.5) cannot have four real zeros at most -1,
# and an antiderivative of the smoothed law f * f_3 can.
import numpy as np, mpmath as mp, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import szego_criterion as S
plt.rcParams.update({"font.family": "serif", "font.size": 9, "mathtext.fontset": "cm"})
mp.mp.dps = 60; f = S.poly_from_gam([1.0, 10.0, 10.5]); h = S.conv(f, S.fpow(3)); h = [c/h[0] for c in h]
def anti(p): return np.array([0.0]+[float(c)/(i+1) for i, c in enumerate(p)])
def val(P, x): return sum(c*x**i for i, c in enumerate(P))
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.7))
for a, p, name in ((ax[0], f, r"antiderivative of $f$"), (ax[1], h, r"antiderivative of $f*f_3$")):
    P = anti(p); r = sorted(float(mp.re(z)) for z in S.roots(p)); x = np.linspace(r[0]-0.25*(r[-1]-r[0])-0.5, -0.6, 2000); crit = [val(P, t) for t in r]
    lo = max(crit[0], crit[2]); hi = min(crit[1], val(P, -1.0)); c = -(lo+hi)/2 if lo <= hi else -lo; y = val(P, x)+c
    a.plot(x, y, "k-", lw=1); a.axhline(0, color="0.6", lw=0.6); a.axvline(-1, color="0.6", lw=0.6, ls=":"); a.plot(r, [val(P, t)+c for t in r], "ko", ms=3, mfc="w")
    s = max(abs(val(P, t)+c) for t in r+[-1.0]); a.set_ylim(-2.5*s, 2.5*s); a.set_xlabel(r"$x=\cos\theta$"); a.set_title(name+(": all zeros at most $-1$" if lo <= hi else ": a zero right of $-1$ or two complex"), fontsize=8)
plt.tight_layout(); plt.savefig("../figures/szego_smoothing.pdf")
