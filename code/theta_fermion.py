# Theta kernels and fermion determinants on SU(2).
# (a) det(1 + alpha U) = 1 + 2 alpha cos(theta) + alpha^2: laws with real zeros are products of fermion determinants.
# (b) Howe duality (dual Cauchy identity): the character coefficients of prod_m det(1 + alpha_m U) are sums of Schur functions of alpha,
#     a_n = sum_q e_q (e_{q+n-1} - e_{q+n+1}), with e_q the elementary symmetric functions.
# (c) The matching expansion of a convolution: f * g = sum over matchings M of the d + D flavours of
#     2^|M|/(|M|+1)! * prod_unmatched (1 + a_i^2) * prod_pairs a_i a_j * x^(number of pairs joining a flavour of f to one of g).
# (d) The ellipse formula: (f*g)(x) is the average of f(s) g(t) over the ellipse s^2 + t^2 - 2 x s t <= 1 - x^2.
import numpy as np, itertools, math
from numpy.polynomial import polynomial as P
def law(al):
    p = np.array([1.0])
    for a in al: p = P.polymul(p, [1+a*a, 2*a])
    return p
def chebU(n): 
    U = [np.array([1.0]), np.array([0.0, 2.0])]
    for k in range(2, n): U.append(P.polysub(P.polymul([0.0, 2.0], U[-1]), U[-2]))
    return U[:n]
def to_char(p):
    d = len(p)-1; U = chebU(d+1); r = np.array(p, float); b = np.zeros(d+1)
    for n in range(d, -1, -1): b[n] = r[n]/U[n][n]; r[:n+1] -= b[n]*U[n]
    return b
def conv(p, q):
    a, b = to_char(p), to_char(q); m = min(len(a), len(b)); U = chebU(m); return sum(a[n]*b[n]/(n+1)*np.pad(U[n], (0, m-len(U[n]))) for n in range(m))
def elem(al):
    e = np.array([1.0])
    for a in al: e = P.polymul(e, [1.0, a])
    return e
def matchings(items):
    if not items: yield []; return
    first, rest = items[0], items[1:]
    for m in matchings(rest): yield m
    for i, other in enumerate(rest):
        for m in matchings(rest[:i]+rest[i+1:]): yield [(first, other)]+m
def matching_expansion(al, be):
    a = list(al)+list(be); d = len(al); n = len(a); out = np.zeros(min(len(al), len(be))+1)
    for M in matchings(list(range(n))):
        used = {i for pr in M for i in pr}; w = 2.0**len(M)/math.factorial(len(M)+1)
        for i in range(n):
            if i not in used: w *= 1+a[i]**2
        cross = 0
        for i, k in M: w *= a[i]*a[k]; cross += (i < d) != (k < d)
        out[cross] += w
    return out
if __name__ == "__main__":
    rng = np.random.default_rng(4); al = rng.uniform(0.1, 1, 3); be = rng.uniform(0.1, 1, 4)
    th = 0.7; print("(a) det(1 + alpha U) against 1 + 2 alpha cos + alpha^2: %.2e" % abs(np.linalg.det(np.eye(2)+al[0]*np.diag([np.exp(1j*th), np.exp(-1j*th)]))-(1+2*al[0]*np.cos(th)+al[0]**2)))
    a = to_char(law(al)); e = np.pad(elem(al), (0, 8)); schur = [sum(e[q]*(e[q+n-1]-e[q+n+1]) for q in range(len(al)+1)) for n in range(1, len(al)+2)]
    print("(b) character coefficients:", np.round(a, 6), " from Schur functions:", np.round(schur, 6))
    c1 = conv(law(al), law(be)); c2 = matching_expansion(al, be); print("(c) convolution by characters:", np.round(c1, 6)); print("    by the matching expansion:", np.round(c2, 6), " (%d matchings)" % sum(1 for _ in matchings(list(range(len(al)+len(be))))))
    x = -0.3; f, g = law(al), law(be); N = 1200; s = (np.arange(N)+0.5)*2/N-1; S, T = np.meshgrid(s, s); inside = S*S+T*T-2*x*S*T <= 1-x*x
    print("(d) ellipse average at x = %.1f: %.5f   convolution: %.5f" % (x, (P.polyval(S, f)*P.polyval(T, g))[inside].mean(), P.polyval(x, c1)))
    z = 0.37+0.0j; lam = al[0]*be[0]; t = np.linspace(0, 400, 2000001); wgt = (1+t)**-4.0; I = [np.trapezoid(wgt*t**k, t) for k in range(3)]
    print("(e) one flavour each: U(2) average 1 + ab x + a^2 b^2, coefficients", np.round([1, lam, lam**2], 5), " from the flavour integral", np.round([1, 2*lam*I[1]/I[0], lam**2*I[2]/I[0]], 5))
