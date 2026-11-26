# The coalescence kernel is hypergeometric: Q(j | a, b) = C(a, j) C(b+2, b-j) / C(a+b+2, b), the law of the size of the intersection
# of two uniformly random subsets of sizes a and b of a set of size a + b + 2.
# 1. Exact check against the character identity u_n(a+1) u_n(b+1) = sum_j Q(j | a, b) u_n(j+1), for a, b <= 25 and all n.
# 2. The drift law of the exact ladder (k = 4, fourfold convolution) at large counts:
#       mean' = mean - 13/16 + (15/16) r,    r' = 5/16 + r/16,    r = variance / mean of the count,
#    against one exact rung from point masses (r = 0), heat kernels (r = 1/3) and Wilson laws (r = 1).
import numpy as np, sys
from fractions import Fraction as F
from math import comb, factorial as fa
sys.path.insert(0, "code")
def u(n, J): return F(fa(J-1)*fa(J+1), fa(J-n)*fa(J+n)) if n <= J else F(0)
def Q(a, b): return [F(comb(a, j)*comb(b+2, b-j), comb(a+b+2, b)) for j in range(min(a, b)+1)]
if __name__ == "__main__":
    ok = all(sum(q*u(n, j+1) for j, q in enumerate(Q(a, b))) == u(n, a+1)*u(n, b+1) for a in range(26) for b in range(a, 26) for n in range(1, a+2))
    print("hypergeometric kernel reproduces the convolution of characters exactly for all a <= b <= 25:", ok)
    import count_lyapunov as L
    from scipy.special import gammaln
    from scipy.optimize import brentq
    j = L.j
    def stats(pi): m = (pi*j).sum(); return m, ((pi*j*j).sum()-m*m)/m
    def after(pi):                        # count law after one rung, by two hypergeometric coalescences of the bond moved law
        nu = np.where(pi > 0, pi, 0)*np.exp(-L.logZ); nu /= nu.sum(); out = nu.copy()
        for _ in range(3): out = np.convolve(out, nu)[:L.MM+1]
        lp = np.where(out > 0, np.log(np.where(out > 0, out, 1))+L.logZ, -np.inf); p = np.exp(lp-lp.max()); p /= p.sum(); uu = (L.U@p)**4
        return 2*(L.n*uu).sum()-2
    cases = {}
    for a in (20, 60): cases["point mass at %d" % a] = L.pad([0]*a+[1])
    for beta in (4.0, 8.0, 16.0):
        t = 2/beta; N = int(60/t)+60; C = np.cosh((2*np.arange(1, N+1)-1)*t/4); f = lambda g: np.sum(1/(C-g)); nu = np.array([1.0])
        for m in range(N-40): p = 1/brentq(f, C[m]*(1+1e-13)+1e-13, C[m+1]*(1-1e-13)-1e-13); nu = np.convolve(nu, [1-p, p])[:L.MM+1]
        nu = np.pad(nu, (0, L.MM+1-nu.size)); cases["heat kernel, beta = %g" % beta] = L.pad(nu*np.exp(L.logZ-L.logZ.max()))
    for k in (8.0, 32.0):
        l = j*np.log(k)-gammaln(j+1)+L.logZ; cases["Wilson, kappa = %g" % k] = L.pad(np.exp(l-l.max()))
    for name, pi in cases.items():
        m, r = stats(pi); print("%-24s mean %7.3f  r = var/mean %.3f  drift predicted %+.4f  exact %+.4f" % (name, m, r, -13/16+15*r/16, after(pi)-m))
