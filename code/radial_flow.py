# The radial flow for convolutions on SU(2).  For f, g with real zeros at most -1 and f_lam(s) = f(lam s):
#   Phi_lam = lam^2 (f_lam * g_lam) = int_0^lam 2 r (f_r o g_r) dr,   with o the convolution on the circle (real zeros by Malo's theorem),
#   f o g = f*g + ((Ef)*g + f*(Eg))/2 with E = x d/dx,
#   and as lam -> 0 the zeros scale like 1/lam^2, with limit polynomial sum_k f_k g_k (y/2)^k/(k+1) (real zeros: Schur, then Laguerre).
#   and (1 + (d+D)/2) f*g = f o g + (1/2) sum_i (f without its i th factor)*g + (1/2) sum_k f*(g without its k th factor).
# Sign criterion: if each of these d + D lower convolutions alternates in sign at the zeros of f o g, then f*g has only real zeros,
# interlaced by those of f o g.
# Checked here: the identities; the sign condition; that the zeros of f_lam * g_lam are real and strictly interlaced by those of f_lam o g_lam along a grid
# of lam; and that they move monotonically.
import mpmath as mp, random
import szego_criterion as S
mp.mp.dps = 80
def chebT(n):
    T = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]]
    for k in range(2, n+1): a = [mp.mpf(0)]+[2*t for t in T[-1]]; b = T[-2]+[mp.mpf(0)]*(len(a)-len(T[-2])); T.append([x-y for x, y in zip(a, b)])
    return T[:n+1]
def to_T(p):
    d = len(p)-1; T = chebT(d); r = p[:]; t = [mp.mpf(0)]*(d+1)
    for n in range(d, -1, -1):
        t[n] = r[n]/T[n][n]
        for i, c in enumerate(T[n]): r[i] -= t[n]*c
    return t
def circ(p, q):
    a, b = to_T(p), to_T(q); m = min(len(a), len(b)); T = chebT(m-1); out = [mp.mpf(0)]*m
    for n in range(m):
        w = a[n]*b[n]*(1 if n == 0 else mp.mpf(1)/2)
        for i, c in enumerate(T[n]): out[i] += w*c
    return out
dil = lambda p, lam: [c*lam**i for i, c in enumerate(p)]
euler = lambda p: [i*c for i, c in enumerate(p)]
def add(*ps): n = max(len(p) for p in ps); return [sum(p[i] if i < len(p) else 0 for p in ps) for i in range(n)]
def real_roots(p):
    r = S.roots(p); assert all(abs(mp.im(z)) < mp.mpf(10)**-20*(1+abs(z)) for z in r), "complex zero"; return sorted(mp.re(z) for z in r)
if __name__ == "__main__":
    random.seed(7); f = S.poly_from_gam([1.3, 2.0, 7.0]); g = S.poly_from_gam([1.1, 1.6, 4.0, 30.0])
    lhs = circ(f, g); rhs = add(S.conv(f, g), [c/2 for c in S.conv(euler(f), g)], [c/2 for c in S.conv(f, euler(g))]); print("identity f o g = f*g + (Ef*g + f*Eg)/2: largest difference %.1e" % float(max(abs(x-y) for x, y in zip(lhs, rhs))))
    X = [mp.mpf(k+mp.mpf(1)/2)/200 for k in range(200)]; quad = [sum(2*r*circ(dil(f, r), dil(g, r))[i] for r in X)/200 for i in range(4)]; print("flow identity, f*g against the radial integral (midpoint rule): %s  %s" % ([mp.nstr(c, 8) for c in S.conv(f, g)], [mp.nstr(c, 8) for c in quad]))
    Q = [f[k]*g[k]/(mp.mpf(2)**k*(k+1)) for k in range(4)]; lam = mp.mpf(1)/200; print("flat limit: zeros of Q %s; lam^2 times zeros of f_lam*g_lam at lam = 1/200 %s" % ([mp.nstr(z, 7) for z in real_roots(Q)], [mp.nstr(lam**2*z, 7) for z in real_roots(S.conv(dil(f, lam), dil(g, lam)))]))
    pairs = cases = 0; worst = mp.inf; mono = True
    for trial in range(60):
        d = random.randint(2, 6); D = random.randint(d, 8); mode = trial % 3
        gam = lambda n: [1+10**random.uniform(-3, 2) for _ in range(n)] if mode < 2 else [random.choice([1.0, 1.02, 3.0, 50.0])]*n
        f, g = S.poly_from_gam(gam(d)), S.poly_from_gam(gam(D)); prev = None; pairs += 1
        for k in range(1, 41):
            lam = mp.mpf(k)/40; zs = real_roots(S.conv(dil(f, lam), dil(g, lam))); zc = real_roots(circ(dil(f, lam), dil(g, lam))); cases += 1
            assert all(zs[i] < zc[i] and (i+1 == len(zs) or zc[i] < zs[i+1]) for i in range(len(zs))), "interlacing fails"
            worst = min(worst, min(min((zc[i]-zs[i])/abs(zs[i]) for i in range(len(zs))), min([(zs[i+1]-zc[i])/abs(zs[i+1]) for i in range(len(zs)-1)] or [mp.inf])))
            if prev is not None: mono = mono and all(a < b for a, b in zip(prev, zs))
            prev = zs
    print("%d pairs, %d values of lam: zeros real and strictly interlaced in every case; smallest relative gap %.2e; zeros increasing in lam: %s" % (pairs, cases, float(worst), mono))
    random.seed(11); tot = ok = 0; err = mp.mpf(0)
    for trial in range(120):
        d = random.randint(2, 6); D = random.randint(2, 7); mode = trial % 3
        gam = lambda n: [1+10**random.uniform(-3, 2) for _ in range(n)] if mode < 2 else [random.choice([1.0, 1.02, 3.0, 50.0])+random.uniform(0, 1e-3) for _ in range(n)]
        gf, gg = gam(d), gam(D); f, g = S.poly_from_gam(gf), S.poly_from_gam(gg); Psi = circ(f, g); psi = sorted(real_roots(Psi), reverse=True)
        members = [S.conv(S.poly_from_gam(gf[:i]+gf[i+1:]), g) for i in range(d)]+[S.conv(f, S.poly_from_gam(gg[:k]+gg[k+1:])) for k in range(D)]
        lhs = [(1+mp.mpf(d+D)/2)*c for c in S.conv(f, g)]; rhs = add(Psi, *[[c/2 for c in M] for M in members]); err = max(err, max(abs(x-y) for x, y in zip(lhs, rhs))/max(abs(x) for x in lhs))
        for M in members: tot += 1; ok += all(mp.sign(S.pv(M, x)) == (-1)**k for k, x in enumerate(psi))
    print("decomposition identity: largest relative error %.1e; sign condition holds for %d of %d lower convolutions" % (float(err), ok, tot))
