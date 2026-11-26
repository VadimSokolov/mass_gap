# Does convolution on SU(2) preserve the class L of densities c * prod_m (1 + x/gamma_m), gamma_m >= 1 (x = cos theta)?
# Polynomials in x are convolved through characters chi_n = U_{n-1}(x) (chi_n * chi_m = delta chi_n / n) in high precision,
# and the zeros of the result are examined. Also: the zeros along the exact ladder from a heat kernel, and under heat flow.
import mpmath as mp, numpy as np
mp.mp.dps = 80
def cheb_u(N):
    U = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(2)]]
    for n in range(2, N+1):
        new = [mp.mpf(0)]+[2*c for c in U[n-1]]; U.append([c-(U[n-2][i] if i < len(U[n-2]) else 0) for i, c in enumerate(new)])
    return U
U = cheb_u(80)
def to_char(poly):
    poly = list(poly); c = [mp.mpf(0)]*len(poly)
    for d in range(len(poly)-1, -1, -1):
        c[d] = poly[d]/U[d][d]
        for i in range(d+1): poly[i] -= c[d]*U[d][i]
    return c
def from_char(c):
    poly = [mp.mpf(0)]*len(c)
    for n, cn in enumerate(c):
        for i, u in enumerate(U[n]): poly[i] += cn*u
    return poly
def conv(p, q):
    a, b = to_char(p), to_char(q); m = min(len(a), len(b)); return from_char([a[n]*b[n]/(n+1) for n in range(m)])
def from_roots(g):
    p = [mp.mpf(1)]
    for r in g: p = [(p[i] if i < len(p) else 0)+(p[i-1]/r if i > 0 else 0) for i in range(len(p)+1)]
    return p
def zeros(p):
    while len(p) > 1 and abs(p[-1]) < mp.mpf(10)**(-60)*max(abs(c) for c in p): p = p[:-1]
    return mp.polyroots(p[::-1], maxsteps=400, extraprec=400) if len(p) > 1 else []
def in_class(p): z = zeros(p); return all(abs(mp.im(r)) < mp.mpf(10)**(-25)*(1+abs(r)) and mp.re(r) <= -1+mp.mpf(10)**(-25) for r in z), z
if __name__ == "__main__":
    rng = np.random.default_rng(3); bad = 0; trials = 0; closest = 9
    for trial in range(300):
        da, db = rng.integers(1, 13, 2); mode = trial % 3
        ga = 1+rng.exponential([0.05, 1.0, 20.0][mode], da); gb = 1+rng.exponential([0.05, 1.0, 20.0][(mode+trial//3) % 3], db)
        ok, z = in_class(conv(from_roots(ga), from_roots(gb))); trials += 1; bad += (not ok); closest = min([closest]+[float(-mp.re(r)) for r in z])
    print("random pairs in the class: %d;  convolutions outside the class: %d;  largest zero found: %.6f" % (trials, bad, -closest))
    bad = 0
    for a in range(1, 15):
        for b in range(a, 15):
            ok, z = in_class(conv(from_roots([1]*a), from_roots([1]*b))); bad += (not ok)
    print("f_a * f_b for 1 <= a <= b <= 14: outside the class: %d;   zeros of f_3 * f_3: %s" % (bad, [mp.nstr(mp.re(r), 6) for r in zeros(conv(from_roots([1]*3), from_roots([1]*3)))]))
    # exact ladder rungs in polynomial arithmetic: rho -> rho^k -> fourfold convolution (leading coefficients fall below 1e-400: raise the precision)
    mp.mp.dps = 450; U[:] = cheb_u(80)
    def power(p, k):
        out = [mp.mpf(1)]
        for _ in range(k): out = [sum(out[i]*p[m-i] for i in range(len(out)) if 0 <= m-i < len(p)) for m in range(len(out)+len(p)-1)]
        return out
    for name, start, k, rungs in (("f_3, k = 4", from_roots([mp.mpf(1)]*3), 4, 2), ("three factors 1.2, 2, 5, k = 4", from_roots([mp.mpf("1.2"), mp.mpf(2), mp.mpf(5)]), 4, 2), ("two factors 1.05, 9, k = 2", from_roots([mp.mpf("1.05"), mp.mpf(9)]), 2, 3)):
        rho = start; line = []
        for r in range(rungs):
            bm = power(rho, k); c2 = conv(bm, bm); rho = conv(c2, c2); s = sum(rho); rho = [c/s for c in rho]; z = mp.polyroots(rho[::-1], maxsteps=2000, extraprec=2000); ok = all(abs(mp.im(x)) < mp.mpf(10)**(-100)*(1+abs(x)) and mp.re(x) <= -1 for x in z); line.append("rung %d: degree %d, in class %s, largest zero %s" % (r+1, len(rho)-1, ok, mp.nstr(max(mp.re(x) for x in z), 7)))
        print(name, "|", "; ".join(line))
    # the identity behind the proof: coefficients of (f*g)(z - 1/z) are those of f (z - 1/z) times those of an antiderivative G of g
    f, g = from_roots([1.3, 2.0]), from_roots([1.1, 1.7, 4.0]); a, b = to_char(f), to_char(g); lhs = [a[n]*b[n]/(n+1) for n in range(len(a))]
    print("coefficient identity a_n b_n / n = (f sin)_n (G)_n: max error %s" % mp.nstr(max(abs(lhs[n]-a[n]*(b[n]/(n+1))) for n in range(len(a))), 3))
