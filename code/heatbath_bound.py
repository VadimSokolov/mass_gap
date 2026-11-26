# Rigorous upper bounds for the 2 x 2 Wilson loop in SU(2) by integrating four loop links that share no plaquette:
#   W22 <= < prod_l u(beta |S_l|) > <= u(2(d-1) beta)^4,   u(x) = I_2(x)/I_1(x),  S_l the staple sum of link l.
# Usage: beta D L
import sys, numpy as np, su2_core as C
from scipy.special import ive
u = lambda x: ive(2, x)/ive(1, x)
beta, D, L = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); C.rng = np.random.default_rng(5); U = C.cold(L, D)
for i in range(100): U = C.sweep(U, beta)
w11, w22, mh = [], [], []
for i in range(40):
    for j in range(4): U = C.sweep(U, beta)
    rho = [u(beta*np.linalg.norm(C.staple(U, mu), axis=-1)) for mu in range(D)]
    for mu in range(D):
        for nu in range(D):
            if mu == nu: continue
            r = rho[mu]; bottom = r*np.roll(r, -1, axis=mu); mh.append((bottom*np.roll(bottom, -2, axis=nu)).mean())
            if mu < nu: w11.append(C.loop(U, mu, nu, 1, 1)[..., 0].mean()); w22.append(C.loop(U, mu, nu, 2, 2)[..., 0].mean())
k = 2**(D-2)
print("d=%d beta %5.1f  W11^4 %.4f  W22 %.4f  <prod rho> %.4f  u(%d beta)^4 %.4f   [bare ladder u(%d beta)^4 %.4f]" % (D, beta, np.mean(w11)**4, np.mean(w22), np.mean(mh), 2*(D-1), u(2*(D-1)*beta)**4, k, u(k*beta)**4), flush=True)
