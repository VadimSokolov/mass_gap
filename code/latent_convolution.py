# Is the family of class functions E[(1+x)^M] on SU(2) closed under convolution? Expand (1+x)^a * (1+x)^b (convolution on the group) in the
# basis (1+x)^j, in exact rational arithmetic, and test the signs of the coefficients. Characters: chi_n = U_{n-1}(x), chi_n * chi_m = delta chi_n / n.
from fractions import Fraction as F
import numpy as np
def cheb_u(N):                      # U_0..U_N as coefficient lists in x
    U = [[F(1)], [F(0), F(2)]]
    for n in range(2, N+1):
        new = [F(0)]+[2*c for c in U[n-1]]; U.append([c-(U[n-2][i] if i < len(U[n-2]) else 0) for i, c in enumerate(new)])
    return U
def to_char(poly, U):               # polynomial in x -> coefficients on chi_1..chi_{deg+1}
    poly = list(poly); c = [F(0)]*len(poly)
    for d in range(len(poly)-1, -1, -1): c[d] = poly[d]/U[d][d]; poly = [p-c[d]*(U[d][i] if i < len(U[d]) else 0) for i, p in enumerate(poly)]
    return c
def binom_pow(a):                   # (1+x)^a
    p = [F(1)]
    for _ in range(a): p = [(p[i] if i < len(p) else 0)+(p[i-1] if i > 0 else 0) for i in range(len(p)+1)]
    return p
def to_shift(poly):                 # polynomial in x -> coefficients in (1+x)^j, by substituting x = z - 1
    out = [F(0)]*len(poly)
    for d, c in enumerate(poly):
        for i in range(d+1): out[i] += c*F((-1)**(d-i))*F(__import__("math").comb(d, i))
    return out
N = 30; U = cheb_u(N); bad = 0; worst = None
for a in range(0, N+1):
    for b in range(a, N+1):
        ca, cb = to_char(binom_pow(a), U), to_char(binom_pow(b), U); m = min(len(ca), len(cb)); conv = [ca[n]*cb[n]/(n+1) for n in range(m)]
        poly = [F(0)]*m
        for n in range(m):
            for i, u in enumerate(U[n]): poly[i] += conv[n]*u
        sh = to_shift(poly)
        if min(sh) < 0: bad += 1; worst = worst or (a, b, [float(s) for s in sh])
print("pairs tested:", (N+1)*(N+2)//2, " pairs with a negative coefficient:", bad)
if worst: print("first failure", worst)
ca = to_char(binom_pow(2), U); print("check a = b = 2:", [str(s) for s in to_shift([sum(ca[n]*ca[n]/(n+1)*(U[n][i] if i < len(U[n]) else 0) for n in range(3)) for i in range(3)])])
# The coalescence kernel: f_a = (1+x)^a / Z_a are probability densities with respect to Haar measure, and f_a * f_b = sum_j Q(j | a, b) f_j.
from math import comb
def Z(a): return F(comb(2*a+2, a+1), 2**a*(a+2))          # integral of (1+x)^a against Haar: Catalan(a+1) / 2^a
def Q(a, b):
    ca, cb = to_char(binom_pow(a), U), to_char(binom_pow(b), U); m = min(len(ca), len(cb)); poly = [F(0)]*m
    for n in range(m):
        for i, u in enumerate(U[n]): poly[i] += ca[n]*cb[n]/(n+1)*u
    sh = to_shift(poly); return [sh[j]*Z(j)/(Z(a)*Z(b)) for j in range(m)]
print("Z(0..3) =", [str(Z(a)) for a in range(4)])
for a, b in ((1, 1), (2, 2), (4, 4), (8, 8), (16, 16), (8, 16), (4, 16)):
    q = Q(a, b); mean = sum(j*float(x) for j, x in enumerate(q)); print("a = %2d b = %2d  sum Q = %s  mean j = %.3f  harmonic ab/(a+b) = %.3f  Q = %s" % (a, b, sum(q), mean, a*b/(a+b), np.round([float(x) for x in q[:6]], 3)))
ok = all(sum(j*x for j, x in enumerate(Q(a, b))) == F(a*b, a+b+2) for a in range(1, 21) for b in range(a, 21)); print("mean of j equals ab/(a+b+2) exactly for all 1 <= a <= b <= 20:", ok)
