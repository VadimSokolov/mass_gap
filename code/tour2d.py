"""Figure tour2d.pdf and the table of Chapter ch:tour (ch07b_tour.tex): two dimensional SU(2), Wilson plaquette law e^{beta cos theta} on Haar.
Formulas: u_n = I_n(beta)/I_1(beta); mean plaquette u_2; gap ratio = -(2 beta/3) log u_2 (lattice gap -(g^2 l beta/4) log u_2 over the continuum gap 3 g^2 l/8);
bound = arcsinh(x)/x with x = 3/(2 beta). Deterministic, no random numbers. Run from code/: python3 tour2d.py"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ive
u = lambda n, b: ive(n, b)/ive(1, b)
ratio = lambda b: -(2*b/3)*np.log(u(2, b)); bound = lambda b: np.arcsinh(1.5/b)/(1.5/b)
print("beta   u2     u3     -log u2  gap ratio  bound")
for b in (0.5, 1, 2, 4, 8, 16, 32):
    print("%g & %.3f & %.3f & %.3f & %.3f & %.3f \\\\" % (b, u(2, b), u(3, b), -np.log(u(2, b)), ratio(b), bound(b)))
bs = np.logspace(np.log10(0.2), np.log10(64), 400)
plt.rcParams.update({"font.size": 10, "axes.spines.top": True})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
for i, n in enumerate(range(2, 8)): ax[0].plot(bs, u(n, bs), color=plt.cm.viridis(i/6), label="$u_{%d}$" % n)
ax[0].set_xscale("log"); ax[0].set_xlabel(r"$\beta$"); ax[0].set_ylabel("character coefficients"); ax[0].legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
ax[1].plot(bs, ratio(bs), color="#c0392b", label="exact gap ratio"); ax[1].plot(bs, bound(bs), "--", color="#1f4e79", label=r"bound $\mathrm{arcsinh}(x)/x$")
ax[1].axhline(1, color="k", lw=0.6); ax[1].set_xscale("log"); ax[1].set_xlabel(r"$\beta$"); ax[1].set_ylabel("gap in continuum units"); ax[1].legend(frameon=False, fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig("../figures/tour2d.pdf"); print("peak of gap ratio: %.4f at beta %.2f" % (ratio(bs).max(), bs[ratio(bs).argmax()]))
