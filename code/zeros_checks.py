"""Checks for chapter 'zeros' (c10_zeros.tex), previously without a script:
 (1) the Gegenbauer identity 2 I_1(r)/r = sum_n (4n/(kappa lambda)) I_n(kappa) I_n(lambda) chi_n(x), r^2 = kappa^2+lambda^2+2 kappa lambda x ("twelve figures"),
 (2) the symbol ((zeta+x)^n * g)(x) of the Borcea Branden criterion: exact for g = 1+x, n = 2, and the random trial count ("720 trials, 559"),
 (3) the Laguerre interpolation of the character coefficients b_m/m of g: complex zeros in how many of 300 random g ("235 of 300").
Convolution on SU(2): f = sum a_n chi_n, g = sum b_n chi_n gives f*g = sum a_n b_n / n chi_n.  The protocols of (2) and (3) are fixed here (seeds 1 and 2);
the chapter does not state them, so the counts are those of THIS protocol."""
import numpy as np
from fractions import Fraction as Fr
from math import comb, factorial
from mpmath import mp, mpf, besseli, sin, acos, sqrt as msqrt
mp.dps = 40
# ---- (1) -------------------------------------------------------------------------------------------------------------
def gegenbauer(kappa, lam, x, N=60):
    th = acos(x); r = msqrt(kappa**2+lam**2+2*kappa*lam*x)
    lhs = 2*besseli(1, r)/r
    rhs = sum(4*n/(kappa*lam)*besseli(n, kappa)*besseli(n, lam)*sin(n*th)/sin(th) for n in range(1, N+1))
    return lhs, rhs
# ---- (2) -------------------------------------------------------------------------------------------------------------
def x_pow_chi(k):                       # x^k in the chi basis, list index m-1 for chi_m;  x chi_m = (chi_{m+1}+chi_{m-1})/2
    v = {1: Fr(1)}
    for _ in range(k):
        w = {}
        for m, c in v.items():
            w[m+1] = w.get(m+1, 0)+c/2
            if m > 1: w[m-1] = w.get(m-1, 0)+c/2
        v = w
    return v
def chi_polys(M):                       # U_0..U_{M-1}
    U = [[Fr(1)], [Fr(0), Fr(2)]]
    for k in range(2, M):
        t = [Fr(0)]+[2*c for c in U[-1]]; p = U[-2]+[Fr(0)]*(len(t)-len(U[-2])); U.append([a-b for a, b in zip(t, p)])
    return U
def g_chi(gam):                         # chi coefficients of g = prod (1 + x/gamma_i), exact
    p = [Fr(1)]
    for c in gam: p = [a+b for a, b in zip(p+[Fr(0)], [Fr(0)]+[y/c for y in p])]
    tot = {}
    for k, ck in enumerate(p):
        for m, c in x_pow_chi(k).items(): tot[m] = tot.get(m, 0)+ck*c
    return tot
def symbol(gam, n):
    """Coefficients (ascending in zeta) of ((zeta+x)^n * g)(x), each a polynomial in x (list ascending), exact."""
    b = g_chi(gam); U = chi_polys(max(n, max(b))+3); out = []
    for j in range(n+1):                # coefficient of zeta^j: comb(n,j) x^(n-j)
        v = x_pow_chi(n-j); poly = [Fr(0)]*(len(U[-1]))
        for m, a in v.items():
            if m in b:
                for i, cc in enumerate(U[m-1]): poly[i] += comb(n, j)*a*b[m]/m*cc
        out.append(poly)
    return out
def run_symbol(trials=720, seed=1):
    rng = np.random.default_rng(seed); bad = 0
    for _ in range(trials):
        d = int(rng.integers(1, 7)); n = int(rng.integers(2, 7)); gam = [Fr(1)+Fr(float(10**rng.uniform(-2, 3))).limit_denominator(10**6) for _ in range(d)]
        S = symbol(gam, n); x = complex(rng.uniform(-3, 3), rng.uniform(0.1, 3))   # x in the open upper half plane
        co = [sum(float(c)*x**i for i, c in enumerate(pj)) for pj in S]
        roots = np.roots(co[::-1])
        if (roots.imag >= 0).any(): bad += 1
    return bad
# ---- (3) -------------------------------------------------------------------------------------------------------------
def complex_zero_interp(gam):
    """G' = g with G = sum_J gamma_J (1+x)^J, G(-1) = 0.  c_m(G) = sum_J gamma_J 2^-J (2J)!/(Gamma(J+1+m)Gamma(J+1-m)), and b_m/m = 2 c_m (up to a constant).
    B(s) = sum_J gamma_J 2^-J (2J)! /(Gamma(J+1+s)Gamma(J+1-s)) = P(s^2)/(Gamma(Jm+1+s)Gamma(Jm+1-s)) with P a polynomial in u = s^2 whose roots are the zeros of B
    other than the real ones +-(Jm+1+k).  B has only real zeros iff every root u of P is real and >= 0.  Returns True if B has a non real zero."""
    p = [Fr(1)]
    for c in gam: p = [a+b for a, b in zip(p+[Fr(0)], [Fr(0)]+[y/c for y in p])]
    # g(x) = sum_k p_k x^k, x = (1+x) - 1, so g = sum_J q_J (1+x)^J; G = sum_J q_J (1+x)^(J+1)/(J+1)
    Jg = len(p)-1; q = [Fr(0)]*(Jg+1)
    for k, pk in enumerate(p):
        for J in range(k+1): q[J] += pk*comb(k, J)*(-1)**(k-J)
    gamJ = {J+1: q[J]/(J+1) for J in range(Jg+1)}; Jm = max(gamJ)
    P = np.polynomial.Polynomial([0.0])
    for J, gj in gamJ.items():
        term = np.polynomial.Polynomial([float(gj/2**J*factorial(2*J))])
        for k in range(J+1, Jm+1): term = term*np.polynomial.Polynomial([float(k*k), -1.0])   # (k+s)(k-s) = k^2 - u
        P = P+term
    u = P.roots()
    return bool(((np.abs(u.imag) > 1e-7*np.maximum(1, np.abs(u))) | (u.real < -1e-9)).any())
def run_interp(trials=300, seed=2):
    rng = np.random.default_rng(seed); cnt = 0
    for _ in range(trials):
        d = int(rng.integers(2, 7)); gam = [Fr(1)+Fr(float(10**rng.uniform(-2, 3))).limit_denominator(10**6) for _ in range(d)]
        cnt += complex_zero_interp(gam)
    return cnt
if __name__ == "__main__":
    for kap, lam, x in [(mpf("1.3"), mpf("0.7"), mpf("0.3")), (mpf(3), mpf(5), mpf("-0.6")), (mpf("0.2"), mpf("9"), mpf("0.95"))]:
        l, r = gegenbauer(kap, lam, x)
        print("Gegenbauer kappa=%s lambda=%s x=%s: lhs %s rhs %s relative difference %.2e" % (kap, lam, x, mp.nstr(l, 20), mp.nstr(r, 20), float(abs(l-r)/abs(l))))
    S = symbol([Fr(1)], 2)
    print("symbol for g = 1+x, n = 2, exact (coefficients of zeta^0, zeta^1, zeta^2, each ascending in x):", [[str(c) for c in pj if True] for pj in S])
    print("random symbol trials with a root zeta in the closed upper half plane (x in the upper half plane): %d of 720" % run_symbol())
    print("Laguerre interpolation of b_m/m has a non real zero: %d of 300" % run_interp())
