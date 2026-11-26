# Decrement per doubling of the lattice spacing, a(b) + c(b)/beta, for the ladder with general scale factor b.
import numpy as np, sys
sys.path.insert(0, "code")
from mk_verify import make_R
def coeffs(b, M=6000, NMAX=1500):
    R = make_R(M, NMAX, b=b); bb = np.array([20., 30, 45, 70, 100, 150])
    d = np.array([x-R(x) for x in bb])/np.log2(b)
    A = np.column_stack([np.ones_like(bb), 1/bb, 1/bb**2]); a, c, e = np.linalg.lstsq(A, d, rcond=None)[0]
    return a, c
if __name__ == "__main__":
    out = []
    for b in [1.2, 1.4, 1.6, 1.8, 1.925, 2.0, 2.2, 2.5, 3.0]:
        a, c = coeffs(b); out.append((b, a, c, (1-b**-2)/(3*np.log2(b))))
        print("b=%.3f  a=%.5f (theory %.5f)  c=%.5f" % (b, a, out[-1][3], c), flush=True)
    np.save("/tmp/ladder_b.npy", np.array(out))
