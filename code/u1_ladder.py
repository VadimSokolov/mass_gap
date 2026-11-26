# The abelian control: the Migdal Kadanoff ladder for compact U(1) in four dimensions, heat kernel (Villain) laws.
# Law on the circle: f_t(theta) = sum_n exp(-n^2 t/2) exp(i n theta). Bond moving: f^k. Decimation: coefficients to the power b^2.
# Coordinate beta = 1/t, read from the first coefficient u_1 = exp(-t/2). Compare with SU(2), where the decrement tends to 1/4.
import numpy as np, sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "code"); from mk_verify import make_R
from mpmath import mp, mpf, exp, cos, quad, log, pi, jtheta
mp.dps = 60
def R_u1(beta, k=4, D=4):
    t = 1/mpf(beta); q = exp(-t/2)
    f = lambda th: jtheta(3, th/2, q)                       # sum_n q^(n^2) e^{i n theta}
    num = quad(lambda th: f(th)**k*cos(th), [0, pi/2, pi]); den = quad(lambda th: f(th)**k, [0, pi/2, pi])
    return 1/(-2*D*log(num/den))
bs = [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0]; d1 = [mpf(b)-R_u1(b) for b in bs]
for b, d in zip(bs, d1): print("U(1) beta %.1f  decrement %.3e   exp(-pi^2 beta) = %.3e" % (b, float(d), float(exp(-pi**2*b))))
R = make_R(3000, 600); b2 = np.geomspace(0.2, 8, 40); d2 = [b-R(b) for b in b2]
fig, ax = plt.subplots(figsize=(5.6, 3.4))
ax.semilogy(b2, d2, color="#c0392b", label="SU(2): tends to $1/4$"); ax.semilogy(bs, [float(abs(d)) for d in d1], "o-", color="#1f4e79", label="U(1): exponentially small")
ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"decrement $\beta-R(\beta)$"); ax.set_ylim(1e-18, 2); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("figures/u1_ladder.pdf")
