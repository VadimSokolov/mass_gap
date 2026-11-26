# The count E[1/(1 - tr U/2)] - 2 of the plaquette and of the 2 x 2 loop on an 8^4 lattice, SU(2), Wilson action.
# The estimator integrates one link exactly: given the rest, the link is W V^-1 with W of class density exp(kappa tr W/2), kappa = beta*|staples|,
# and E[1/(1 - tr(W Y)/2)] = 2 sum_n (I_n(kappa)/I_1(kappa)) chi_n(Y); this has finite variance, which the plain average has not.
# Compared with the counts of a Wilson law and of a heat kernel having the same fundamental coefficient. Usage: beta [configs]
import sys, numpy as np
from scipy.special import ive
from scipy.optimize import brentq
import su2_core as C
N = np.arange(1, 90)
def est(kap, y0):
    al = np.arccos(np.clip(y0, -1+1e-12, 1-1e-12)); out = np.zeros_like(kap); i1 = ive(1, kap)
    for n in N: out += ive(n, kap)/i1*np.sin(n*al)/np.sin(al)
    return 2*out
def wilson_count(w): k = brentq(lambda k: ive(2, k)/ive(1, k)-w, 1e-6, 400); return 2*(N*ive(N, k)/ive(1, k)).sum()-2
def heat_count(w): t = -4*np.log(w)/3; return 2*(N*np.exp(-(N*N-1)*t/4)).sum()-2
def measure(U, beta):
    D = 4; c1, c2, w1, w2 = [], [], [], []
    for mu in range(D):
        A = C.staple(U, mu); k = np.linalg.norm(A, axis=-1); Vc = C.cj(A/k[..., None]); Um = U[..., mu, :]
        for nu in range(D):
            if nu == mu: continue
            for R in (1, 2):
                P = C.loop(U, mu, nu, R, R); Y = C.qm(Vc, C.qm(C.cj(Um), P)); e = est(beta*k, Y[..., 0]); (c1 if R == 1 else c2).append(e.mean()-2); (w1 if R == 1 else w2).append(P[..., 0].mean())
    return np.mean(c1), np.mean(c2), np.mean(w1), np.mean(w2)
if __name__ == "__main__":
    beta = float(sys.argv[1]); nc = int(sys.argv[2]) if len(sys.argv) > 2 else 30; C.rng = np.random.default_rng(int(100*beta)); U = C.cold(8, 4)
    for i in range(80): U = C.sweep(U, beta)
    out = []
    for i in range(nc):
        for j in range(4): U = C.sweep(U, beta)
        out.append(measure(U, beta))
    out = np.array(out); m = out.mean(0); se = out.std(0)/np.sqrt(nc)
    print("beta %.2f  W11 %.4f  W22 %.4f | count 1x1 %.4f (%.4f)  Wilson matched %.4f  heat kernel matched %.4f | count 2x2 %.4f (%.4f)  Wilson matched %.4f  heat kernel matched %.4f"
          % (beta, m[2], m[3], m[0], se[0], wilson_count(m[2]), heat_count(m[2]), m[1], se[1], wilson_count(m[3]), heat_count(m[3])))
