# The exact Migdal Kadanoff ladder as a process on integer counts.
# Densities f_j = (1+x)^j / Z_j on SU(2) have character coefficients u_n(J) = (J-1)!(J+1)! / ((J-n)!(J+n)!), J = j+1, 1 <= n <= J.
# Coalescence theorem: u_n(A) u_n(B) = sum_J Q_J u_n(J) with Q >= 0 and mean of j equal to ab/(a+b+2).
# Part 1 (exact arithmetic): compute Q by the triangular solve, test positivity, the mean, and pairwise total positivity of Q(j | a, b).
# Part 2 (floating point): run the exact ladder (k = 4, b^2 = 4) on count laws from a heat kernel and check monotonicity along the orbit.
import numpy as np, sys, pickle, os
from fractions import Fraction as F
from math import factorial as fa, comb
def u(n, J): return F(fa(J-1)*fa(J+1), fa(J-n)*fa(J+n)) if n <= J else F(0)
def Qrow(a, b):
    A, B = min(a, b)+1, max(a, b)+1; rhs = [u(n, A)*u(n, B) for n in range(1, A+1)]; q = [F(0)]*(A+1)
    for J in range(A, 0, -1):                                         # u_n(J) = 0 for n > J: solve from the top
        q[J] = (rhs[J-1]-sum(q[K]*u(J, K) for K in range(J+1, A+1)))/u(J, J)
    return q[1:]                                                      # index j = J - 1
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
if __name__ == "__main__":
    Q = {(a, b): Qrow(a, b) for a in range(NMAX+1) for b in range(a, NMAX+1)}; get = lambda j, a, b: (Q[(min(a, b), max(a, b))][j] if j <= min(a, b) else F(0))
    print("pairs a <= b <= %d: %d;  negative entries: %d;  mean = ab/(a+b+2) in all: %s" % (NMAX, len(Q), sum(x < 0 for q in Q.values() for x in q), all(sum(j*x for j, x in enumerate(q)) == F(a*b, a+b+2) for (a, b), q in Q.items())))
    bad_ja = sum(get(j, a, b)*get(j+1, a+1, b) < get(j+1, a, b)*get(j, a+1, b) for b in range(NMAX+1) for a in range(NMAX) for j in range(NMAX))
    bad_ab = sum(get(j, a, b)*get(j, a+1, b+1) < get(j, a+1, b)*get(j, a, b+1) for j in range(NMAX) for a in range(NMAX) for b in range(NMAX))
    st = sum(sum(get(i, a, b) for i in range(j+1)) < sum(get(i, a+1, b) for i in range(j+1)) for b in range(NMAX+1) for a in range(NMAX) for j in range(NMAX))
    print("violations of TP2 in (j, a) at fixed b: %d;  of TP2 in (a, b) at fixed j: %d;  of stochastic monotonicity in a: %d" % (bad_ja, bad_ab, st))
    print("Q(. | 3, 3) =", [str(x) for x in Q[(3, 3)]], "  Q(. | 1, b) = [3/(b+3), b/(b+3)]:", all(Q[(1, b)] == [F(3, b+3), F(b, b+3)] for b in range(1, NMAX+1)))
    # ---- Part 2: the exact ladder on counts, k = 4 and b^2 = 4, from a heat kernel
    import mpmath as mp
    from scipy.optimize import brentq
    from scipy.special import gammaln
    mp.mp.dps = 60; MM = 110
    lu = lambda n, J: gammaln(J)+gammaln(J+2)-gammaln(J-n+1)-gammaln(J+n+1)
    Uf = np.array([[np.exp(lu(n, J)) if n <= J else 0.0 for J in range(1, MM+2)] for n in range(1, MM+2)])       # U[n-1, J-1]
    Um = mp.matrix(MM+1, MM+1)
    for n in range(1, MM+2):
        for J in range(n, MM+2): Um[n-1, J-1] = mp.factorial(J-1)*mp.factorial(J+1)/(mp.factorial(J-n)*mp.factorial(J+n))
    logZ = np.array([gammaln(2*m+3)-2*gammaln(m+2)-np.log(m+2)-m*np.log(2) for m in range(MM+1)])
    def decimate(pi):                                                # pi -> law of the fourfold convolution, through characters
        un = mp.matrix([mp.mpf(float(x))**4 for x in Uf@pi]); q = [mp.mpf(0)]*(MM+1)
        for J in range(MM+1, 0, -1): q[J-1] = (un[J-1]-sum(q[K-1]*Um[J-1, K-1] for K in range(J+1, MM+2)))/Um[J-1, J-1]
        return np.array([float(x) for x in q])
    def bond(pi, k=4):
        nu = pi*np.exp(-logZ); nu /= nu.sum(); out = nu.copy()
        for _ in range(k-1): out = np.convolve(out, nu)[:MM+1]
        p = out*np.exp(logZ); return p/p.sum(), nu
    for beta in (1.0, 2.0, 3.0, 4.0):
        t = 2/beta; N = int(60/t)+50; C = np.cosh((2*np.arange(1, N+1)-1)*t/4); f = lambda g: np.sum(1/(C-g))
        p = np.array([1/brentq(f, C[m]*(1+1e-13)+1e-13, C[m+1]*(1-1e-13)-1e-13) for m in range(N-40)]); nu = np.array([1.0])
        for q in p: nu = np.convolve(nu, [1-q, q])[:MM+1]
        nu = np.pad(nu, (0, MM+1-nu.size)); pi = nu*np.exp(logZ); pi /= pi.sum(); print("beta %.1f" % beta); orbit = []
        for rung in range(8):
            j = np.arange(MM+1); u2 = (pi*j/(j+3)).sum(); bm, nu0 = bond(pi); new = decimate(bm); s = pi > 1e-12; lc = np.all(nu0[1:-1][(nu0[2:] > 1e-13)]**2 >= (nu0[:-2]*nu0[2:])[(nu0[2:] > 1e-13)]*(1-1e-9))
            both = (pi > 1e-10) & (new > 1e-10); r = (new/np.where(pi > 0, pi, 1))[both]; lr = np.all(np.diff(r) <= 1e-9*r[:-1])
            print("  rung %d  u_2 %.6f  mean count %.3f  min pi %.1e  nu_0 log concave %s  next law below in likelihood ratio order %s  mass beyond 100: %.1e" % (rung, u2, (pi*j).sum(), new.min(), lc, lr, bm[100:].sum()))
            orbit.append(pi[:16].copy()); pi = np.clip(new, 0, None); pi /= pi.sum()
        np.save("code/count_orbit_%d.npy" % beta, np.array(orbit))
