# Lowest order lattice perturbation theory for Wilson loops:
# -log W(R,R) = (g^2 C_F / 2) (1/V) sum_k F_R(k_mu) F_R(k_nu) (khat_mu^2 + khat_nu^2)/khat^2,  F_R(k) = sin^2(kR/2)/sin^2(k/2).
# Returns w(R)/w(1), the weak coupling limit of log W_RR / log W_11, on an L^d lattice (zero mode dropped).
import numpy as np
def ratios(d, L, Rs=(1, 2, 4)):
    k = 2*np.pi*np.arange(L)/L; grids = np.meshgrid(*([k]*d), indexing="ij")
    kh2 = [4*np.sin(g/2)**2 for g in grids]; tot = sum(kh2); tot[(0,)*d] = np.inf
    out = []
    for R in Rs:
        def F(g):
            with np.errstate(divide="ignore", invalid="ignore"):
                f = np.sin(g*R/2)**2/np.sin(g/2)**2
            return np.where(np.isfinite(f), f, R*R)
        out.append((F(grids[0])*F(grids[1])*(kh2[0]+kh2[1])/tot).mean())
    out = np.array(out); return out/out[0], out[0]
if __name__ == "__main__":
    for d, Ls in [(4, (8, 16, 32)), (3, (16, 32, 64, 128))]:
        for L in Ls:
            r, w1 = ratios(d, L); print("d=%d L=%3d  w11 coefficient %.5f (2/d = %.5f)  w22/w11 %.4f  w44/w11 %.4f" % (d, L, w1, 2/d, r[1], r[2]))
