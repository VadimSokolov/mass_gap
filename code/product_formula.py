# Product formula for the SU(2) heat kernel in c = cos(theta), and monotonicity of its roots in t.
# K_t = (e^{t/4}/2) C(q) P'(c),  P(c) = prod_n (1 + a_n^2 + 2 a_n c),  a_n = q^(2n-1),  q = e^{-t/4}.
# Zeros of P: c_n = -cosh((2n-1)t/4).  Zeros of P' interlace: K_t = A(t) prod_m (1 + c/gamma_m(t)).
import numpy as np
from mpmath import mp, mpf, cosh, findroot, exp, sin, cos, nsum, inf
mp.dps = 40
def gammas(t, M=8, NT=60):
    cn = [-cosh((2*n-1)*mpf(t)/4) for n in range(1, NT+1)]
    S = lambda x: sum(1/(x-c) for c in cn)
    out = []
    for m in range(M):
        lo, hi = cn[m+1], cn[m]; d = (hi-lo)*mpf(10)**(-25)
        out.append(-findroot(S, (lo+d, hi-d), solver="anderson", tol=1e-30))
    return out
def K(t, th):
    return nsum(lambda n: n*exp(-(n*n-1)*mpf(t)/4)*sin(n*th)/sin(th), [1, inf])
if __name__ == "__main__":
    for t in [0.5, 1, 2, 4]:
        g = gammas(t, M=40, NT=120)
        prod = lambda th: np.prod([1+cos(th)/x for x in g])
        r = [K(t, th)/prod(th) for th in [0.3, 1.0, 2.0, 2.9]]
        print("t=%.1f gamma_1..3 = %s  ratio spread %.2e" % (t, [float(x) for x in g[:3]], float(max(r)/min(r)-1)))
    ts = np.linspace(0.2, 6, 30); G = np.array([[float(x) for x in gammas(t, 5)] for t in ts])
    print("all gamma_m increasing in t:", bool((np.diff(G, axis=0) > 0).all()))
    np.save("/tmp/gammas.npy", np.column_stack([ts, G]))
